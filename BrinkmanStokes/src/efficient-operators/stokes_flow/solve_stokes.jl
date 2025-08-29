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
	println("Mesh size: $(msize) ($(2^(msize+1))×$(2^(msize+1)) grid)")
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
	println("  Original system dimensions:")
	println("    A matrix: $(size(A))")
	println("    B matrix: $(size(B))")
	println("    f vector: $(length(f))")
	println("    g vector: $(length(g))")
	println("    Boundary nodes: $(length(bound))")

	@time begin
		Ast, Bst, fst, gst = flowbc(A, B, f, g, xy, bound, domain)
	end

	println("  Constrained system dimensions:")
	println("    Ast matrix: $(size(Ast))")
	println("    Bst matrix: $(size(Bst))")
	println("    fst vector: $(length(fst))")
	println("    gst vector: $(length(gst))")

	np = length(gst)
	rhs = [fst; gst]

	println("  Constrained velocity DOFs: $(size(Ast, 1))")
	println("  Pressure DOFs: $(np)")

	# Assemble saddle-point system matrix
	println("Assembling saddle-point system...")
	@time begin
		K = [Ast Bst'; Bst spzeros(np, np)]
	end

	println("  System matrix: $(size(K, 1))×$(size(K, 2))")
	println("  Nonzeros: $(nnz(K))")
	println("  Sparsity: $(round(nnz(K)/length(K)*100, digits=2))%")
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
	tol = 1e-6                       # Convergence tolerance
	maxIter = 100                    # Maximum iterations

	# Add timeout for large problems
	timeout_seconds = 300  # 5 minutes timeout
	if size(K, 1) > 10000
		timeout_seconds = 1800  # 30 minutes for very large problems
		println("  Large system detected - extending timeout to $(timeout_seconds/60) minutes")
	end

	println("GMRES Configuration:")
	println("  Restart: $(restrt)")
	println("  Tolerance: $(tol)")
	println("  Max iterations: $(maxIter)")
	println("  Preconditioning: $(use_mg ? "Multigrid" : "None")")
	println()

	# Analyze system properties before solving
	println("Analyzing system matrix properties...")
	@time begin
		println("  Matrix dimensions: $(size(K))")
		println("  Matrix nnz: $(nnz(K))")
		println("  Matrix density: $(round(nnz(K)/length(K)*100, digits=4))%")

		# Check matrix symmetry
		K_symmetric = ishermitian(K)
		println("  Matrix is symmetric: $(K_symmetric)")

		# Sample some matrix-vector products for timing
		test_vec = randn(size(K, 2))
		@time K_test = K * test_vec
		println("  Sample matrix-vector product completed")

		# Check for obvious issues
		rhs_norm = norm(rhs)
		println("  RHS norm: $(rhs_norm)")
		if rhs_norm < 1e-14
			println("  WARNING: RHS norm is very small - may indicate singular system")
		end

		# Check diagonal entries
		K_diag = diag(K)
		min_diag = minimum(abs.(K_diag[K_diag.!=0]))
		max_diag = maximum(abs.(K_diag))
		println("  Diagonal entries - min: $(min_diag), max: $(max_diag)")
		if min_diag / max_diag < 1e-12
			println("  WARNING: Large diagonal ratio suggests ill-conditioning")
		end

		# Memory usage estimate
		matrix_memory_mb = (nnz(K) * 16 + size(K, 1) * 8) / 1024^2  # 16 bytes per nonzero + 8 per index
		println("  Estimated matrix memory: $(round(matrix_memory_mb, digits=2)) MB")
	end
	println()

	# Create custom GMRES callback for detailed logging
	iteration_count = 0
	last_log_time = time()

	function gmres_callback(x, r)
		global iteration_count, last_log_time
		iteration_count += 1
		current_time = time()

		if iteration_count == 1 || iteration_count % 10 == 0 || (current_time - last_log_time) > 30.0
			residual_norm = norm(r)
			println("  GMRES iteration $(iteration_count): residual = $(residual_norm)")
			last_log_time = current_time
		end

		# Check for stagnation
		if iteration_count > 50 && length(resvec) > 10
			recent_progress = resvec[end-9] / resvec[end]
			if recent_progress < 1.01  # Less than 1% improvement over 10 iterations
				println("  WARNING: GMRES appears to be stagnating (progress ratio: $(recent_progress))")
			end
		end

		return false  # Continue iteration
	end

	# Solve the linear system with timeout
	println("Solving linear system...")
	println("Starting GMRES iteration...")
	println("Timeout set to: $(timeout_seconds) seconds")
	start_solve_time = time()

	# Wrap GMRES in a timeout mechanism
	solve_task = @async begin
		gmres(K, rhs, restrt; tol = tol, maxIter = maxIter, M = M, out = 1)
	end

	# Check for timeout
	xst, flag, err, iter, resvec = nothing, -99, Inf, 0, Float64[]
	solve_completed = false

	while !solve_completed && (time() - start_solve_time) < timeout_seconds
		if istaskdone(solve_task)
			try
				xst, flag, err, iter, resvec = fetch(solve_task)
				solve_completed = true
			catch e
				println("ERROR: GMRES failed with exception: $e")
				break
			end
		else
			sleep(1.0)  # Check every second
			elapsed = time() - start_solve_time
			if elapsed > 30 && elapsed % 30 < 1  # Log every 30 seconds after first 30 seconds
				println("  GMRES still running... elapsed time: $(round(elapsed, digits=1))s")
			end
		end
	end

	if !solve_completed
		println("ERROR: GMRES timed out after $(timeout_seconds) seconds")
		println("Consider:")
		println("  1. Using multigrid preconditioning (set use_mg=true)")
		println("  2. Reducing mesh size")
		println("  3. Checking matrix conditioning")
		error("GMRES solver timeout - system may be ill-conditioned")
	end

	solve_time = time() - start_solve_time
	println("Total solve time: $(round(solve_time, digits=2)) seconds")

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
