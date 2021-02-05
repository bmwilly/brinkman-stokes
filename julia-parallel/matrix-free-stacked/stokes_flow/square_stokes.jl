# load required packages and functions
include("stokes_q2p1.jl")
include("../helpers/helper_functions.jl")
include("../grids/q2p1grid.jl")
include("../grids/cavity_domain.jl")

###SQUARE_STOKES set up flow problem in unit square domain
function square_stokes(msize)

    # generate Q2 grid for square cavity
    cavity_grid = cavity_domain(msize)
    grid = q2p1grid(cavity_grid)

    # stokes q2-p1 matrix generator
    mv = cavity_grid["mv"]
    stokes_grid = merge(grid, {"mv" => mv})
    stokes_mats = stokes_q2p1(stokes_grid)

    bounds = {
      "bound" => cavity_grid["bound"],
      "bndxy" => cavity_grid["bndxy"],
      "bnde" => cavity_grid["bnde"],
      "obs" => cavity_grid["obs"]
    }

    order = 2; dim = 2;
    nelems = [2^(msize)]
    m = Mesh.Hexmesh(tuple(repeat(nelems, 1, dim)...), Xform.identity)
    dof = prod([m.nelems...] * order + 1)
    Mesh.set_order(m, order);
    refel = Refel(m.dim, order);
    dof = prod([m.nelems...] * order + 1);
    ne = prod([m.nelems...]);
    # storage for indices and values
    NP = (order + 1)^m.dim;
    NPNP = NP * NP;
    bdy = Mesh.get_boundary_node_indices(m, order);

    params = {
      "mesh" => m,
      "order" => order,
      "dof" => dof,
      "ne" => ne,
      "NP" => NP,
      "bdy" => bdy,
      "refel" => refel,
      "mv" => mv
    }

    kparams = merge(stokes_mats, stokes_grid, bounds, params, {"msize" => msize})

end
