###Q2P1GRID Q2-P1 element grid generator
# input
#   grid          grid for a domain with 5 properties
#     x             x coordinate vector
#     y             y coordinate vector
#     xy            nodal coordinate vector
#     mv            Q2 macroelement mapping matrix
#     bound         boundary vertex vector
# output
#   grid_out
#     xyp           centroid coordinate vector
#     ee            element edge connection matrix
function q2p1grid(grid)

    x = grid["x"]; y = grid["y"]; xy = grid["xy"]
    mv = grid["mv"]; bound = grid["bound"]

    ## centroid coordinate vector
    xx = xy[:, 1]
    yy = xy[:, 2]
    nvtx = length(xx)
    nel = length(mv[:, 1])

    ## recompute mid-side points in the case of stretched grids
    # y-direction
    yv = yy
    ny = length(y)

    for k in 2:2:ny
        yold = y[k]
        ynew = 0.5 * (y[k + 1] + y[k - 1])
        l = findall(yy == yold)
        yv[l] .= ynew
        y[k] = ynew
    end

    # x-direction
    xv = xx
    nx = length(x)

    for k in 2:2:nx
        xold = x[k]
        xnew = 0.5 * (x[k + 1] + x[k - 1])
        l = findall(xx == xold)
        xv[l] .= xnew
        x[k] = xnew
    end

    xy = [xv yv]

    # centroid coordinates
    xc = zeros(nel, 1)
    yc = zeros(nel, 1)
    for ielem in 1:nel
        xc[ielem] = mean(xx[mv[ielem, 1:4]])
        yc[ielem] = mean(yy[mv[ielem, 1:4]])
    end

    xyp = [xc yc]

    # compute edge to edge connection array ee
    np = nel
    # initialize global matrices
    adj = spzeros(nvtx, nvtx)
    ee = []

    return grid_out = {
        "x" => x,
        "y" => y,
        "xy" => xy,
        "xyp" => xyp,
        "ee" => ee,
        "mv" => mv,
    }

end
