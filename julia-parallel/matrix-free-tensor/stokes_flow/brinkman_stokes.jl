using LinearOperators
using Distances
reload("grids/q2p1grid.jl")
reload("stokes_flow/stokes_q2p1.jl")
reload("grids/channel_domain.jl")
reload("stokes_flow/ho_afun.jl")
reload("stokes_flow/ho_gfun.jl")
include("../helpers/helper_functions.jl")
include("../../julia-homg/Basis.jl")
include("../../julia-homg/Hexmesh.jl")
include("../../julia-homg/Xform.jl")
include("../../julia-homg/Grids.jl")
include("../../julia-homg/Tensor.jl")
include("../../julia-homg/Refel.jl")

function brinkman_stokes(msize)

  channel_grid = channel_domain(msize) # Q2 grid for channel domain
  stokes_grid = q2p1grid(channel_grid)
  stokes_mats = stokes_q2p1(stokes_grid) # stokes element matrices

  bounds = {
    "bound" => channel_grid["bound"],
    "bndxy" => channel_grid["bndxy"],
    "bnde" => channel_grid["bnde"],
    "obs" => channel_grid["obs"]
  }

  # brinkman obstacles
  centers = [
    0.1   0.5;
    0.15  0.7;
    0.18  0.1;
    0.2   0.8;
    0.45  0.55;
    0.5   0.35;
    0.65  0.65;
    0.75  0.4;
    0.77  0.33;
    0.8   0.42;
  ]

  msize = channel_grid["msize"]; order = 2; dim = 2;
  nelems = [2^(msize-1)]
  m = Mesh.Hexmesh(tuple(repmat(nelems, 1, dim)...), Xform.identity)
  dof = prod([m.nelems...]*order + 1)
  Mesh.set_order(m,order);
  refel = Refel( m.dim, order );
  dof = prod([m.nelems...]*order + 1);
  ne = prod([m.nelems...]);
  # storage for indices and values
  NP = (order+1)^m.dim;
  NPNP = NP * NP;
  bdy = Mesh.get_boundary_node_indices(m, order);
  pts =  Mesh.element_nodes(m, 1, refel);
  (detJac, Jac) = Mesh.geometric_factors(m, refel, pts);
  brinkman_pts = Mesh.brinkman_tensor(pts, centers)
  eMat_stiffness = Mesh.element_stiffness_brinkman(m, 1, refel, detJac, Jac, brinkman_pts);
  eMat_mass = Mesh.element_mass_brinkman(m, 1, refel, detJac, brinkman_pts);

  params = {
    "mesh" => m,
    "order" => order,
    "dof" => dof,
    "ne" => ne,
    "NP" => NP,
    "bdy" => bdy
  }
  params_stiffness = merge(params, {"eMat" => eMat_stiffness})
  params_mass = merge(params, {"eMat" => eMat_mass})

  A = LinearOperator(2dof, Float64, u -> ho_afun(u, params_stiffness))
  G = LinearOperator(2dof, Float64, u -> ho_gfun(u, params_mass))

  mats = merge(stokes_mats, stokes_grid, bounds, params_stiffness, {"A" => A, "G" => G})

end
