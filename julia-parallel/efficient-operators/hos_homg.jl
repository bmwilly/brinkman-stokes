# reload("helpers/helper_functions.jl")
# reload("../julia-homg/Basis.jl")
# reload("../julia-homg/Hexmesh.jl")
# reload("../julia-homg/Xform.jl")
# reload("../julia-homg/Grids.jl")
# reload("../julia-homg/Tensor.jl")
# reload("../julia-homg/Refel.jl")
include("helpers/helper_functions.jl")
include("../julia-homg/Basis.jl")
include("../julia-homg/Hexmesh.jl")
include("../julia-homg/Xform.jl")
include("../julia-homg/Grids.jl")
include("../julia-homg/Tensor.jl")
include("../julia-homg/Refel.jl")

function hos_homg(order, msize, dim)
  # dim = 2
  nelems = [2^msize]

  m = Mesh.Hexmesh(tuple(repmat(nelems, 1, dim)...), Xform.identity)
  dof = prod([m.nelems...]*order + 1)
  K,M,iK = Mesh.assemble_poisson(m, order)
  k1,k2 = size(K)
  A = [K spzeros(k1,k2); spzeros(k1,k2) K]
  tic()
  for cnt = 1:100; u = rand(2dof); w = A*u; end
  etoc = toc()
end
