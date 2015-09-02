@everywhere using ParallelSparseMatMul
@everywhere include("helpers/helper_functions.jl")
@everywhere include("../julia-homg/Basis.jl")
@everywhere include("../julia-homg/Hexmesh.jl")
@everywhere include("../julia-homg/Xform.jl")
@everywhere include("../julia-homg/Grids.jl")
@everywhere include("../julia-homg/Tensor.jl")
@everywhere include("../julia-homg/Refel.jl")

function hos_homg(order, msize, dim)
  nelems = [2^msize]
  m = Mesh.Hexmesh(tuple(repmat(nelems, 1, dim)...), Xform.identity)
  dof = prod([m.nelems...]*order + 1)
  K,M = Mesh.assemble_poisson(m, order)
  k1,k2 = size(K)
  A = [K spzeros(k1,k2); spzeros(k1,k2) K]
  S = share(A)
  L = operator(A)
  tic()
  for cnt = 1:100
    u = Base.shmem_rand(2dof);
    w = S*u;
  end
  etoc = toc()
end
