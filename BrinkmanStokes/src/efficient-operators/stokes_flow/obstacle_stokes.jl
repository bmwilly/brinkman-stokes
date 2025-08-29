reload("grids/obstacle_domain.jl")
reload("grids/q2p1grid.jl")
reload("stokes_flow/stokes_q2p1.jl")

###OBSTACLE_STOKES sets up flow problem in domain with square obstacle
function obstacle_stokes()

    # generate Q2 grid for obstacle domain
    obstacle_grid = obstacle_domain()
    grid = q2p1grid(obstacle_grid)

    # stokes q2-p1 matrix generator
    stokes_grid = merge(grid, {"mv" => obstacle_grid["mv"]})
    stokes_mats = stokes_q2p1(stokes_grid)

    bounds = {
        "bound" => obstacle_grid["bound"],
        "bndxy" => obstacle_grid["bndxy"],
        "bnde" => obstacle_grid["bnde"],
        "obs" => obstacle_grid["obs"],
    }
    return mats = merge(stokes_mats, grid, bounds)

end
