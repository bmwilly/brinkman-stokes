# @everywhere using ParallelSparseMatMul
@everywhere include("helpers/helper_functions.jl")
# @everywhere include("../julia-homg/Basis.jl")
# @everywhere include("../julia-homg/Hexmesh.jl")
# @everywhere include("../julia-homg/Xform.jl")
# @everywhere include("../julia-homg/Grids.jl")
# @everywhere include("../julia-homg/Tensor.jl")
# @everywhere include("../julia-homg/Refel.jl")
# @everywhere include("stokes_flow/ho_afun.jl")
# using ParallelSparseMatMul
# include("helpers/helper_functions.jl")
# include("../julia-homg/Basis.jl")
# include("../julia-homg/Hexmesh.jl")
# include("../julia-homg/Xform.jl")
# include("../julia-homg/Grids.jl")
# include("../julia-homg/Tensor.jl")
# include("../julia-homg/Refel.jl")
# include("stokes_flow/ho_afun.jl")
# using Basis
# using Mesh
# using Xform
# using Grids
# using Tensor
# using HexMeshGrids

function hos_homg(order, msize, dim)
  # dim = 2
  # @everywhere order = order
  # @everywhere msize = msize
  # @everywhere dim = dim
  # @everywhere nelems = [2^msize]

  nelems = [2^msize]

  m = Mesh.Hexmesh(tuple(repmat(nelems, 1, dim)...), Xform.identity)
  Mesh.set_order(m,order);
  refel = Refel( m.dim, order );
  dof = prod([m.nelems...]*order + 1);
  ne = prod([m.nelems...]);
  # storage for indices and values
  NP = (order+1)^m.dim;
  bdy = Mesh.get_boundary_node_indices(m, order);
  pts =  Mesh.element_nodes(m, 1, refel);
  (detJac, Jac) = Mesh.geometric_factors(m, refel, pts);
  eMat = Mesh.element_stiffness(m, 1, refel, detJac, Jac);

  params = {
    "mesh" => m,
    "order" => order,
    "dof" => dof,
    "ne" => ne,
    "NP" => NP,
    "bdy" => bdy,
    "eMat" => eMat
  }

  # @everywhere A = LinearOperator(2dof, Float64, u -> ho_afun(u, params))

  # tic()
  # for cnt = 1:100
  #   # u = rand(2dof)
  #   u = SharedArray(Float64, 2dof)
  #   w = ho_afun(u, params)
  # end
  # etoc = toc()
end
