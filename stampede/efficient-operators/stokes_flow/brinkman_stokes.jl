# load required packages and functions
using MAT
using Distances
reload("stokes_flow/stokes_brinkman_q2p1.jl")
reload("grids/q2p1grid.jl")
reload("grids/channel_domain.jl")
include("../helpers/helper_functions.jl")
include("../../julia-homg/Basis.jl")
include("../../julia-homg/Hexmesh.jl")
include("../../julia-homg/Xform.jl")
include("../../julia-homg/Grids.jl")
include("../../julia-homg/Tensor.jl")
include("../../julia-homg/Refel.jl")

###SQUARE_STOKES set up flow problem in unit square domain
function brinkman_stokes()

    # generate Q2 grid for square channel
    channel_grid = channel_domain()
    grid = q2p1grid(channel_grid)

    # stokes q2-p1 matrix generator
    mv = channel_grid["mv"]
    stokes_grid = merge(grid, {"mv" => mv})

    xy = stokes_grid["xy"]
    nu = 2length(xy[:, 1])
    K = zeros(nu, nu)
    stokes_mats = stokes_brinkman_q2p1(stokes_grid, K)

    bounds = {
      "bound" => channel_grid["bound"],
      "bndxy" => channel_grid["bndxy"],
      "bnde" => channel_grid["bnde"],
      "obs" => channel_grid["obs"]
    }

    msize = channel_grid["msize"];
    order = int(input("Polynomial order: "))
    # dim = int(input("Dimension: "));
    dim = 2;
    nelems = [2^(msize-1)]
    m = Mesh.Hexmesh(tuple(repmat(nelems, 1, dim)...), Xform.identity)
    dof = prod([m.nelems...]*order + 1)

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
    K,M,iK = Mesh.assemble_poisson_brinkman(m, order, centers)
    k1,k2 = size(K)
    A = [K spzeros(k1,k2); spzeros(k1,k2) K]
    G = [M spzeros(k1,k2); spzeros(k1,k2) M]
    stokes_mats["A"] = A
    stokes_mats["G"] = G

    # keys(mats) =
    # {"A", "B", "G", "Q", "Bx", "By", "f", "g", "x", "y", "xyp", "bound"}
    mats = merge(stokes_mats, grid, bounds)
end