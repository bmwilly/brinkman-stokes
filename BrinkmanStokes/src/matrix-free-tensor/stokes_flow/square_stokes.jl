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
	stokes_mats = stokes_q2p1(grid)

	bounds = Dict(
		"bound" => cavity_grid["bound"],
		"bndxy" => cavity_grid["bndxy"],
		"bnde" => cavity_grid["bnde"],
		"obs" => cavity_grid["obs"],
	)

	kparams = merge(stokes_mats, grid, bounds, Dict("msize" => msize))
	return kparams

end
