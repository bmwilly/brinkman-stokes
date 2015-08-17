using ParallelSparseMatMul
# reload("helpers/helper_functions.jl")
# reload("../julia-homg/Basis.jl")
# reload("../julia-homg/Hexmesh.jl")
# reload("../julia-homg/Xform.jl")
# reload("../julia-homg/Grids.jl")
# reload("../julia-homg/Tensor.jl")
# reload("../julia-homg/Refel.jl")
@everywhere include("helpers/helper_functions.jl")
@everywhere include("../julia-homg/Basis.jl")
@everywhere include("../julia-homg/Hexmesh.jl")
@everywhere include("../julia-homg/Xform.jl")
@everywhere include("../julia-homg/Grids.jl")
@everywhere include("../julia-homg/Tensor.jl")
@everywhere include("../julia-homg/Refel.jl")

function hos_homg(order, msize, dim)
  # dim = 2
  nelems = [2^msize]

  m = Mesh.Hexmesh(tuple(repmat(nelems, 1, dim)...), Xform.identity)
  dof = prod([m.nelems...]*order + 1)
  K,M,iK = Mesh.assemble_poisson(m, order)
  k1,k2 = size(K)
  A = [K spzeros(k1,k2); spzeros(k1,k2) K]
  S = share(A)
  L = operator(A)
  tic()
  for cnt = 1:100
    u = share(rand(2dof));
    # w = S*u;
    w = L*u;
  end
  etoc = toc()
end
