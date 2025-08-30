# load required packages and functions
using MAT
using Distances

using BrinkmanStokes: Mesh, Xform

include("stokes_brinkman_q2p1.jl")
include("../grids/q2p1grid.jl")
include("../grids/channel_domain.jl")
include("../helpers/helper_functions.jl")


###SQUARE_STOKES set up flow problem in unit square domain
function brinkman_stokes(msize)

	# generate Q2 grid for square channel
	println("Grid generation for channel domain.")
	@time channel_grid = channel_domain(msize)
	@time grid = q2p1grid(channel_grid)

	# stokes q2-p1 matrix generator
	mv = channel_grid["mv"]
	stokes_grid = merge(grid, Dict("mv" => mv))

	xy = stokes_grid["xy"]
	nu = 2length(xy[:, 1])
	# K = zeros(nu, nu)
	# stokes_mats = stokes_brinkman_q2p1(stokes_grid, K)
	stokes_mats = stokes_q2p1(stokes_grid)

	bounds = Dict(
		"bound" => channel_grid["bound"],
		"bndxy" => channel_grid["bndxy"],
		"bnde" => channel_grid["bnde"],
		"obs" => channel_grid["obs"],
	)

	# order = user_input("Polynomial order: ")
	order = 2
	# dim = user_input("Dimension: ");
	dim = 2
	# nelems = [2^(msize-1)]
	nelems = [2^msize]
	m = Mesh.Hexmesh(tuple(repeat(nelems, 1, dim)...), Xform.identity)
	dof = prod([m.nelems...] * order .+ 1)

	# brinkman obstacles
	centers = [
		0.15 0.2
		0.17 0.45
		0.15 0.8
		0.3  0.3
		0.4  0.55
		0.4  0.8
		0.7  0.4
		0.7  0.8
		0.8  0.1
		0.9  0.35
	]
	# centers = [0.33 0.5]

	# spe10
	# x = grid["x"]; y = grid["y"];
	# spe10 = matread("data/spe10.mat")
	# KU = spe10["KU"]; pU = spe10["pU"];
	# # layer = 1;
	# # K = KU[1, 1:((length(x)-1)/2), 1:((length(y)-1)/2), layer];
	# # K = squeeze(K[1,:,:], 1);
	# centers = KU[:] .^(-1)
	# kp = reshape(kp[1:length(x)*length(y)], length(x), length(y))

	# kp = KU[1, 1, :, :]
	# kp = squeeze(kp, 1)
	# kp = squeeze(kp, 1)

	println("Assembling Brinkman matrices...")
	@time K, M, kappa = Mesh.assemble_poisson_brinkman(m, order, centers)
	k1, k2 = size(K)
	println("  Brinkman matrix K: $(k1)×$(k2), nnz: $(nnz(K))")
	@time A = [K spzeros(k1, k2); spzeros(k1, k2) K]
	println("  Block matrix A: $(size(A)), nnz: $(nnz(A))")
	G = [M spzeros(k1, k2); spzeros(k1, k2) M]
	stokes_mats["A"] = A
	stokes_mats["G"] = G

	mats = merge(stokes_mats, grid, bounds, Dict("kappa" => kappa, "msize" => msize))
	return mats
end
