# push!(LOAD_PATH, "/Users/bwilliams/Documents/brinkman-stokes/julia-parallel/julia-homg")
# push!(LOAD_PATH, "/home/bmw313/Documents/brinkman-stokes/julia-parallel/julia-homg")
# @everywhere using Basis
# @everywhere using Mesh
# @everywhere using Xform
# @everywhere using Grids
# @everywhere using Tensor
# @everywhere using HexMeshGrids
# @everywhere require("Refel.jl")

function hos_homg_setup(order, msize, dim)
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
end
