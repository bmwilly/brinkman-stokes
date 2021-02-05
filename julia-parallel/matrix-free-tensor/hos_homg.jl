# reload("helpers/helper_functions.jl")
# reload("../julia-homg/Basis.jl")
# reload("../julia-homg/Hexmesh.jl")
# reload("../julia-homg/Xform.jl")
# reload("../julia-homg/Grids.jl")
# reload("../julia-homg/Tensor.jl")
# reload("../julia-homg/Refel.jl")
using LinearOperators
include("helpers/helper_functions.jl")
include("../julia-homg/Basis.jl")
include("../julia-homg/Hexmesh.jl")
include("../julia-homg/Xform.jl")
include("../julia-homg/Grids.jl")
include("../julia-homg/Tensor.jl")
include("../julia-homg/Refel.jl")
include("stokes_flow/ho_afun.jl")

function hos_homg(order, msize, dim)
  # dim = 2
    nelems = [2^msize]

  # m = Mesh.Hexmesh(tuple(repeat(nelems, 1, dim)...), Xform.identity)
  # dof = prod([m.nelems...]*order + 1)
  # K,M,iK = Mesh.assemble_poisson(m, order)
  # k1,k2 = size(K)
  # A = [K spzeros(k1,k2); spzeros(k1,k2) K]
  # tic()
  # for cnt = 1:100; u = vec(rand(2dof, 1)); w = A*u; end
  # etoc = toc()

    m = Mesh.Hexmesh(tuple(repeat(nelems, 1, dim)...), Xform.identity)
    Mesh.set_order(m, order);
    refel = Refel(m.dim, order);
    dof = prod([m.nelems...] * order + 1);
    ne = prod([m.nelems...]);
  # storage for indices and values
    NP = (order + 1)^m.dim;
    NPNP = NP * NP;
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

    A = LinearOperator(2dof, Float64, u -> ho_afun(u, params))
    tic()
    for cnt = 1:100; u = rand(2dof); w = A * u; end
    etoc = toc()
end
