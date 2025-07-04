push!(LOAD_PATH, "$(homedir())/Documents/brinkman-stokes/julia-parallel/julia-homg")
@everywhere using LinearOperators
@everywhere using ParallelSparseMatMul
@everywhere using Mesh
@everywhere using Xform
@everywhere include("helpers/helper_functions.jl")
@everywhere include("grids/channel_domain.jl")
@everywhere include("../julia-homg/Refel.jl")
@everywhere include("stokes_flow/ho_afun.jl")
@everywhere include("stokes_flow/ho_afun_square.jl")

function hos_homg(order, msize, dim)
	# dim = 2
	nelems = [2^msize]

	# channel_grid = channel_domain(msize) # Q2 grid for channel domain
	# mv = channel_grid["mv"]

	m = Mesh.Hexmesh(tuple(repeat(nelems, 1, dim)...), Xform.identity)
	Mesh.set_order(m, order)
	refel = Refel(m.dim, order)
	dof = prod([m.nelems...] * order .+ 1)
	ne = prod([m.nelems...])
	# storage for indices and values
	NP = (order + 1)^m.dim
	NPNP = NP * NP
	bdy = Mesh.get_boundary_node_indices(m, order)
	pts = Mesh.element_nodes(m, 1, refel)
	(detJac, Jac) = Mesh.geometric_factors(m, refel, pts)
	eMat = Mesh.element_stiffness(m, 1, refel, detJac, Jac)

	idxs = zeros(ne, NP)
	for e ∈ 1:ne
		idx = Mesh.get_node_indices(m, e, order)
		idxs[e, :] = idx
	end

	params = {
		"mesh" => m,
		"order" => order,
		"dof" => dof,
		"ne" => ne,
		"NP" => NP,
		"bdy" => bdy,
		"eMat" => eMat,
		"refel" => refel,
		"mv" => mv,
		"idxs" => idxs,
	}

	# A = LinearOperator(2dof, Float64, u -> ho_afun(u, params))
	# A = LinearOperator(2dof, Float64, u -> ho_afun_square(u, params))
	A = u -> ho_afun_square(u, params)
	# u = Base.shmem_rand(2dof)
	# w = A(u)
	tic()
	for cnt ∈ 1:100
		u = Base.shmem_rand(2dof)
		w = A(u)
	end
	etoc = toc()
end
