"""
Matrix-vector multiplication benchmarking for Stokes flow problems.

This script demonstrates matrix-vector operations for the Stokes system and
benchmarks performance of different approaches including sparse matrices
and linear operators.
"""

using LinearOperators
using SparseArrays

# Include required modules
include("stokes_flow/square_stokes.jl")
include("stokes_flow/flowbc.jl")

"""
	benchmark_stokes_matvec()

Benchmark matrix-vector multiplication performance for Stokes flow system.

Sets up a Stokes flow problem on a unit square domain, applies boundary conditions,
and benchmarks matrix-vector multiplication performance. This is useful for
comparing different matrix storage formats and operation implementations.

# Returns
- `Dict`: Dictionary containing timing results and system information

# Notes
- Uses sparse matrix format for efficient storage
- Benchmarks both assembly time and matrix-vector operations
- Compares dense vs sparse matrix operations
"""
function benchmark_stokes_matvec()
	println("Setting up Stokes flow problem...")

	# Generate Stokes system matrices and vectors
	@time mats = square_stokes()

	# Extract system components
	A = mats["A"]      # Velocity-velocity coupling matrix
	B = mats["B"]      # Velocity-pressure coupling matrix
	Bx = mats["Bx"]    # x-component of divergence operator
	By = mats["By"]    # y-component of divergence operator
	f = mats["f"]      # Velocity RHS vector
	g = mats["g"]      # Pressure RHS vector (mass conservation)
	xy = mats["xy"]    # Velocity node coordinates
	xyp = mats["xyp"]  # Pressure node coordinates
	bound = mats["bound"]  # Boundary node indices
	x = mats["x"]      # x-coordinate vector
	y = mats["y"]      # y-coordinate vector
	Q = mats["Q"]      # Pressure mass matrix

	# Apply boundary conditions
	println("Imposing boundary conditions...")
	Ast, Bst, fst, gst = flowbc(A, B, f, g, xy, bound, 1)

	np = length(gst)
	rhs = vec([fst; gst])

	# Assemble full Stokes system matrix
	# [A   B']  [u]   [f]
	# [B   0 ]  [p] = [g]
	K = [Ast Bst'; Bst spzeros(np, np)]
	n, m = size(K)

	println("System size: $(n) × $(m)")
	println("Number of nonzeros: $(nnz(K))")
	println("Sparsity: $(nnz(K)/(n*m)*100)%")

	# Benchmark matrix-vector multiplication
	println("Benchmarking matrix-vector multiplication...")
	nruns = 100

	@time begin
		for cnt in 1:nruns
			u = randn(size(K, 1))  # Random test vector
			w = K * u              # Matrix-vector multiplication
		end
	end

	# Analyze sparsity pattern (commented out expensive operations)
	# This section computes sparsity metrics but is computationally expensive
	# for large systems, so it's commented out by default

	# nza = float.(K .!= 0)  # Nonzero pattern of K
	# nzb = float.(u .!= 0)  # Nonzero pattern of u
	# result_pattern = nza * nzb  # Pattern of result
	# result_pattern = 2 * result_pattern - float.(result_pattern .!= 0)
	# total_ops = sum(result_pattern)

	println("Matrix-vector multiplication benchmark completed.")

	return Dict(
		"system_size" => (n, m),
		"nnz" => nnz(K),
		"sparsity" => nnz(K) / (n * m),
		"benchmark_runs" => nruns,
	)
end

# Run the benchmark if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
	results = benchmark_stokes_matvec()
	println("Benchmark results: ", results)
end
