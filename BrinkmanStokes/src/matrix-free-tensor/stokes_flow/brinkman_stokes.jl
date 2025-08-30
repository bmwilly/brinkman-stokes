using LinearOperators
using Distances
include("../grids/q2p1grid.jl")
include("../grids/channel_domain.jl")
include("stokes_q2p1.jl")
include("ho_afun.jl")
include("ho_gfun.jl")
include("../helpers/helper_functions.jl")
include("../../julia-homg/Basis.jl")
include("../../julia-homg/Hexmesh.jl")
include("../../julia-homg/Xform.jl")
include("../../julia-homg/Grids.jl")
include("../../julia-homg/Tensor.jl")
include("../../julia-homg/Refel.jl")

function brinkman_stokes(msize)

	channel_grid = channel_domain(msize) # Q2 grid for channel domain
	grid = q2p1grid(channel_grid)
	mv = channel_grid["mv"]
	stokes_grid = merge(grid, Dict("mv" => mv))
	stokes_mats = stokes_q2p1(stokes_grid) # stokes element matrices

	bounds = Dict(
		"bound" => channel_grid["bound"],
		"bndxy" => channel_grid["bndxy"],
		"bnde" => channel_grid["bnde"],
		"obs" => channel_grid["obs"],
	)

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
	# centers = [
	#   0.1   0.5;
	#   0.15  0.7;
	#   0.18  0.1;
	#   0.2   0.8;
	#   0.45  0.55;
	#   0.5   0.35;
	#   0.65  0.65;
	#   0.75  0.4;
	#   0.77  0.33;
	#   0.8   0.42;
	# ]
	# centers = [0.33 0.5]

	order = 2
	dim = 2
	nelems = [2^(msize)]
	m = Mesh.Hexmesh(tuple(repeat(nelems, 1, dim)...), Xform.identity)
	dof = prod([m.nelems...] * order .+ 1)
	Mesh.set_order(m, order)
	refel = Refel(m.dim, order)
	dof = prod([m.nelems...] * order .+ 1)
	ne = prod([m.nelems...])
	# storage for indices and values
	NP = (order + 1)^m.dim
	NPNP = NP * NP
	bdy = Mesh.get_boundary_node_indices(m, order)

	params = Dict(
		"mesh" => m,
		"order" => order,
		"dof" => dof,
		"ne" => ne,
		"NP" => NP,
		"bdy" => bdy,
		"refel" => refel,
		"centers" => centers,
		"mv" => mv,
	)

	A = LinearOperator(2dof, Float64, u -> ho_afun(u, params))
	G = LinearOperator(2dof, Float64, u -> ho_gfun(u, params))

	mats = merge(stokes_mats, stokes_grid, bounds, params, Dict("A" => A, "G" => G))
	return mats

end
