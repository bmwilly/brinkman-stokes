# load required packages and functions
reload("helpers/meshgrid.jl")
reload("stokes_flow/stokes_q2p1.jl")
reload("helpers/input.jl")
reload("grids/q2p1grid.jl")
reload("grids/cavity_domain.jl")

###SQUARE_STOKES set up flow problem in unit square domain
function square_stokes(msize)

    # generate Q2 grid for square cavity
    cavity_grid = cavity_domain(msize)
    grid = q2p1grid(cavity_grid)

    # stokes q2-p1 matrix generator
    stokes_grid = merge(grid, {"mv" => cavity_grid["mv"]})
    stokes_mats = stokes_q2p1(stokes_grid)

    bounds = {
      "bound" => cavity_grid["bound"],
      "bndxy" => cavity_grid["bndxy"],
      "bnde" => cavity_grid["bnde"],
      "obs" => cavity_grid["obs"]
    }

    # keys(mats) =
    # {"A", "B", "G", "Q", "Bx", "By", "f", "g", "x", "y", "xyp", "bound"}
    mats = merge(stokes_mats, grid, bounds)

end
