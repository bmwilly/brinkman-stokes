# load required packages and functions
# include("../helpers/meshgrid.jl")
# include("stokes_q2p1.jl")
# include("../helpers/input.jl")
# require("helpers/meshgrid.jl")
# require("stokes_flow/stokes_q2p1.jl")
reload("stokes_flow/stokes_q2p1.jl")
reload("helpers/helper_functions.jl")
reload("grids/q2p1grid.jl")
reload("grids/cavity_domain.jl")

###SQUARE_STOKES set up flow problem in unit square domain
function square_stokes(msize)

    # generate Q2 grid for square cavity
    cavity_grid = cavity_domain(msize)
    grid = q2p1grid(cavity_grid)

    # stokes q2-p1 matrix generator
    stokes_mats = stokes_q2p1(grid)

    bounds = {
      "bound" => cavity_grid["bound"],
      "bndxy" => cavity_grid["bndxy"],
      "bnde" => cavity_grid["bnde"],
      "obs" => cavity_grid["obs"]
    }

    kparams = merge(stokes_mats, grid, bounds, {"msize" => msize})

end
