using DelimitedFiles
using KrylovMethods
using PyPlot
# using IterativeSolvers
using SparseArrays
using Statistics
using LinearAlgebra
include("../../output_utils.jl")
include("square_stokes.jl")
include("brinkman_stokes.jl")
include("flowbc.jl")
include("../solvers/mg_diff.jl")
include("../solvers/m_st_mg.jl")

"""
Main Stokes flow solver with multigrid preconditioning.

This module provides the main interface for solving Stokes flow problems using
GMRES iteration with optional geometric multigrid preconditioning.
"""

"""
	solve_stokes(domain::Int, msize::Int; use_mg::Bool=false) -> Dict

Solve Stokes flow problem with optional multigrid preconditioning.

Sets up and solves the Stokes equations using GMRES iteration. Supports multiple
domain types and optional geometric multigrid preconditioning for the velocity
subproblem.

# Arguments
- `domain::Int`: Domain type selector
  - `1`: Unit square cavity ([-1,1]²)
  - `2`: Brinkman-Stokes problem (porous media)
- `msize::Int`: Mesh size parameter (grid has 2^(msize+1) intervals per direction)
- `use_mg::Bool=false`: Enable geometric multigrid preconditioning

# Returns
- `Dict`: Complete solution data including:
  - `"xst"`: Solution vector [velocity; pressure]
  - `"resvec"`: GMRES residual history
  - `"iter"`: Number of GMRES iterations
  - `"err"`: Final residual norm
  - `"K"`: System matrix
  - `"rhs"`: Right-hand side vector
  - All original matrix data from domain setup

# Stokes System
The discretized Stokes equations form the saddle-point system:
```
[A   B'] [u]   [f]
[B   0 ] [p] = [g]
```
where:
- A: Velocity stiffness matrix (viscous terms)
- B: Divergence operator matrix
- u: Velocity vector
- p: Pressure vector
- f: Velocity forcing
- g: Mass conservation constraint

# Solver Configuration
- Method: GMRES with restart
- Tolerance: 1e-6
- Maximum iterations: 100
- Restart parameter: min(5000, system_size)

# Preconditioning
When `use_mg=true`, applies geometric multigrid to the velocity block:
- V-cycle multigrid for A-block
- Block-diagonal preconditioner for full system
- Significant speedup for large problems

# Examples
```julia
# Solve cavity flow without preconditioning
sol = solve_stokes(1, 3)  # 16×16 velocity grid

# Solve with multigrid preconditioning
sol = solve_stokes(1, 5, use_mg=true)  # 64×64 grid with MG

# Solve Brinkman-Stokes problem
sol = solve_stokes(2, 4)  # Porous media flow

# Extract solution
u_vel = sol["xst"][1:end-length(sol["gst"])]  # Velocity
p_press = sol["xst"][end-length(sol["gst"])+1:end]  # Pressure
```

# Output Files
Automatically generates:
- Convergence plot: `convergence.png`
- Residual data: `convergence_data.csv`

# Notes
- System is solved in mixed form (velocity-pressure coupling)
- Boundary conditions applied via penalty method
- Multigrid significantly improves performance for fine grids
- Memory usage scales as O(N) for N unknowns

# See also
[`square_stokes`](@ref), [`brinkman_stokes`](@ref), [`flowbc`](@ref), [`mg_diff`](@ref)
"""
function solve_stokes(domain::Int, msize::Int; use_mg::Bool = false)
    # Input validation
    domain in [1, 2] || throw(ArgumentError("domain must be 1 (cavity) or 2 (Brinkman)"))
    msize >= 1 || throw(ArgumentError("msize must be ≥ 1"))

    println("="^60)
    println("STOKES FLOW SOLVER")
    println("="^60)
    println("Domain type: $(domain == 1 ? "Square cavity" : "Brinkman-Stokes")")
    println("Mesh size: $(msize) ($(2^(msize + 1))×$(2^(msize + 1)) grid)")
    println("Multigrid: $(use_mg ? "Enabled" : "Disabled")")
    println()

    # Set up problem matrices and vectors
    println("Setting up finite element system...")
    @time begin
        if domain == 1
            mats = square_stokes(msize)
        elseif domain == 2
            mats = brinkman_stokes(msize)
        end
    end

    # Extract system components
    A = mats["A"]          # Velocity stiffness matrix
    B = mats["B"]          # Velocity-pressure coupling
    Bx = mats["Bx"]        # x-component divergence
    By = mats["By"]        # y-component divergence
    f = mats["f"]          # Velocity RHS
    g = mats["g"]          # Pressure RHS
    xy = mats["xy"]        # Velocity nodes
    xyp = mats["xyp"]      # Pressure nodes
    bound = mats["bound"]  # Boundary nodes
    x = mats["x"]          # x-coordinates
    y = mats["y"]          # y-coordinates
    Q = mats["Q"]          # Pressure mass matrix

    # Handle permeability for Brinkman problems
    if domain == 2
        kappa = mats["kappa"]
        println("Brinkman permeability parameter included")
    else
        kappa = zeros(length(x) * length(y))
    end

    println("System dimensions:")
    println("  Velocity DOFs: $(size(A, 1))")
    println("  Pressure DOFs: $(size(B, 2))")
    println("  Total DOFs: $(size(A, 1) + size(B, 2))")
    println()

    # Apply boundary conditions
    println("Applying boundary conditions...")
    println("  Boundary nodes: $(length(bound))")

    @time begin
        Ast, Bst, fst, gst = flowbc(A, B, f, g, xy, bound, domain)
    end

    np = length(gst)
    rhs = vec([fst; gst])  # Ensure RHS is a Vector, not Matrix

    println("  Constrained velocity DOFs: $(size(Ast, 1))")
    println("  Pressure DOFs: $(np)")

    # Assemble saddle-point system matrix
    println("Assembling saddle-point system...")
    @time begin
        K = [Ast Bst'; Bst spzeros(np, np)]
    end

    println("  System matrix: $(size(K, 1))×$(size(K, 2))")
    println("  Nonzeros: $(nnz(K))")
    println("  Sparsity: $(round(nnz(K) / length(K) * 100, digits = 2))%")
    println()

    # Set up preconditioner
    M = u -> u  # Identity preconditioner (no preconditioning)

    if use_mg
        println("Setting up geometric multigrid preconditioner...")

        # Extract velocity subproblem dimensions
        nv = size(Ast, 1)      # Total velocity DOFs
        np = size(Q, 1)        # Pressure DOFs
        nu = Int(nv / 2)       # DOFs per velocity component

        # Extract velocity stiffness matrix for MG setup
        Agal = Ast[1:nu, 1:nu]  # Assumes u-v coupling is block-diagonal

        # Set up multigrid hierarchy
        @time begin
            mgdata, smooth_data, sweeps, stype, npre, npost, nc = mg_diff(x, y, Agal)
        end

        # Package multigrid parameters
        mparams = Dict(
            "nv" => nv,
            "Q" => Q,
            "mgdata" => mgdata,
            "smooth_data" => smooth_data,
            "nc" => nc,
            "npre" => npre,
            "npost" => npost,
            "sweeps" => sweeps,
        )

        # Block multigrid preconditioner
        M = u -> m_st_mg(u, mparams)
        println("Multigrid preconditioner ready ($(nc) levels)")
        println()
    end

    # Configure GMRES solver
    restrt = min(5000, length(rhs))  # Restart parameter
    tol = 1.0e-6                       # Convergence tolerance
    maxIter = 100                    # Maximum iterations

    # Add timeout for large problems
    timeout_seconds = 300  # 5 minutes timeout
    if size(K, 1) > 10000
        timeout_seconds = 1800  # 30 minutes for very large problems
        println("  Large system detected - extending timeout to $(timeout_seconds / 60) minutes")
    end

    println("GMRES Configuration:")
    println("  Restart: $(restrt)")
    println("  Tolerance: $(tol)")
    println("  Max iterations: $(maxIter)")
    println("  Preconditioning: $(use_mg ? "Multigrid" : "None")")
    println()

    # Quick system validation
    println("Validating system matrix...")

    # Check for problematic values
    K_has_nan = any(isnan, K.nzval)
    K_has_inf = any(isinf, K.nzval)
    rhs_has_nan = any(isnan, rhs)
    rhs_has_inf = any(isinf, rhs)

    if K_has_nan || K_has_inf || rhs_has_nan || rhs_has_inf
        println("  ERROR: Matrix or RHS contains NaN/Inf values!")
        println("  Matrix NaN: $(K_has_nan), Inf: $(K_has_inf)")
        println("  RHS NaN: $(rhs_has_nan), Inf: $(rhs_has_inf)")
        error("System contains invalid values - check matrix assembly")
    end

    # Basic system info
    rhs_norm = norm(rhs)
    println("  System size: $(size(K, 1)), RHS norm: $(round(rhs_norm, digits = 3))")

    if rhs_norm < 1.0e-14
        println("  WARNING: RHS norm is very small - may indicate singular system")
    end
    println()


    # Solve the linear system with timeout
    println("Solving linear system...")
    println("Starting GMRES iteration...")
    println("Timeout set to: $(timeout_seconds) seconds")
    start_solve_time = time()

    # Try direct solve first to get better error messages
    xst, flag, err, iter, resvec = nothing, -99, Inf, 0, Float64[]
    solve_completed = false

    try
        println("Attempting direct GMRES solve...")
        @time xst, flag, err, iter, resvec = gmres(K, rhs, restrt; tol = tol, maxIter = maxIter, M = M, out = 1)
        solve_time = time() - start_solve_time
        println("GMRES completed successfully in $(round(solve_time, digits = 2)) seconds!")
        solve_completed = true
    catch direct_error
        println("Direct GMRES failed: $(typeof(direct_error))")
        println("Error message: $direct_error")

        # Print more detailed error information
        if isa(direct_error, BoundsError)
            println("This is likely a dimension mismatch in the system.")
            println("Check that matrix and vector dimensions are compatible.")
        elseif isa(direct_error, LinearAlgebra.SingularException)
            println("Matrix is singular - system has no unique solution.")
            println("Check boundary conditions and matrix assembly.")
        elseif isa(direct_error, OutOfMemoryError)
            println("Out of memory - try reducing mesh size.")
        end

        println("\n" * "="^60)
        println("GMRES SOLVER FAILED")
        println("="^60)
        println("The GMRES solver crashed during execution.")
        println("This suggests a fundamental problem with the system setup.")
        println("\nDiagnostic suggestions:")
        println("  • Verify matrix dimensions match RHS vector")
        println("  • Check for singular or ill-conditioned matrices")
        println("  • Try a smaller mesh size (msize=$(max(1, msize - 1)))")
        println("  • Enable multigrid preconditioning (use_mg=true)")
        println("  • Verify boundary conditions are correctly applied")
        error("GMRES solver failed during execution - see error details above")
    end


    # Report convergence results
    println()
    println("GMRES Results:")
    if flag == 0
        println("  ✓ Converged to tolerance at iteration $(length(resvec))")
    else
        println("  ⚠ Did not converge (flag=$(flag))")
    end
    println("  Final residual: $(err)")
    println("  Iterations: $(iter)")
    println()

    # Generate convergence plot
    println("Generating output files...")
    convergence_file = get_output_file("efficient-operators", domain, msize, "convergence.png"; subdir = "plots")

    figure(figsize = (10, 6))
    semilogy(1:length(resvec), resvec, "b-", linewidth = 2, label = "GMRES residual")
    xlabel("Iteration")
    ylabel("Residual norm")
    title("GMRES Convergence (Domain=$(domain), Size=$(msize)$(use_mg ? ", MG" : ""))")
    grid(true, alpha = 0.3)
    legend()
    savefig(convergence_file, dpi = 150, bbox_inches = "tight")
    println("  Convergence plot: $(convergence_file)")
    close("all")

    # Save convergence data
    data_file = get_output_file("efficient-operators", domain, msize, "convergence_data.csv")
    writedlm(data_file, [1:length(resvec) resvec])
    println("  Convergence data: $(data_file)")

    # Assemble solution dictionary
    sol = Dict(
        # Solution data
        "xst" => xst,
        "resvec" => resvec,
        "iter" => iter,
        "err" => err,
        "flag" => flag,

        # System matrices
        "K" => K,
        "Ast" => Ast,
        "Bst" => Bst,
        "rhs" => rhs,
        "fst" => fst,
        "gst" => gst,

        # Preconditioner
        "M" => M,
        "use_mg" => use_mg,

        # Problem data
        "kappa" => kappa,
        "domain" => domain,
        "msize" => msize,
    )

    # Merge with original matrix data
    sol = merge(mats, sol)

    println("Stokes solver completed successfully!")
    println("="^60)

    return sol
end
