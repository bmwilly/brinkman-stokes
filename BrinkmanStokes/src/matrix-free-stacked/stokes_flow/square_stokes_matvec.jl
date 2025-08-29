# load required packages and functions
include("../helpers/meshgrid.jl")
include("stokes_q2p1.jl")
include("../helpers/input.jl")

###SQUARE_STOKES set up flow problem in unit square domain
function square_stokes_matvec(msize)

    ## define geometry
    println("Grid generation for cavity domain.")
    n = 2^msize
    np = int(n / 2)
    nel = int(np^2)

    # y-direction
    yy = [(1 / np):(1 / np):1]
    ypos = [0, yy]
    yneg = -yy[length(yy):-1:1]
    y = [yneg, ypos]'
    x = y

    # compute biquadratic element coordinates
    nvtx = (n + 1) * (n + 1)
    (X, Y) = meshgrid(x, y)
    xx = reshape(X', nvtx, 1)
    yy = reshape(Y', nvtx, 1)
    xy = [xx[:] yy[:]]

    kx = 1
    ky = 1
    mel = 0
    mv = zeros(Int64, nel, 9)
    for j in 1:np
        for i in 1:np
            mref = (n + 1) * (ky - 1) + kx
            mel += 1
            nvv = zeros(9)
            nvv[1] = mref
            nvv[2] = mref + 2
            nvv[3] = mref + 2n + 4
            nvv[4] = mref + 2n + 2
            nvv[5] = mref + 1
            nvv[6] = mref + n + 3
            nvv[7] = mref + 2n + 3
            nvv[8] = mref + n + 1
            nvv[9] = mref + n + 2
            mv[mel, 1:9] = nvv[1:9]
            kx += 2
        end
        ky += 2
        kx = 1
    end

    # compute boundary vertices and edges
    # four boundary edges
    k1 = findall(xy[:, 2] .== -1)
    e1 = []
    for k in 1:mel
        if any(mv[k, 5] .== k1)
            e1 = [e1, k]
        end
    end
    ef1 = ones(size(e1))

    k2 = findall((xy[:, 1] .== 1) .& (xy[:, 2] .< 1) .& (xy[:, 2] .> -1))
    e2 = []
    for k in 1:mel
        if any(mv[k, 6] .== k2)
            e2 = [e2, k]
        end
    end
    ef2 = 2 * ones(size(e2))

    k3 = findall(xy[:, 2] .== 1)
    e3 = []
    for k in 1:mel
        if any(mv[k, 7] .== k3)
            e3 = [e3, k]
        end
    end
    ef3 = 3 * ones(size(e3))

    k4 = findall((xy[:, 1] .== -1) .& (xy[:, 2] .< 1) .& (xy[:, 2] .> -1))
    e4 = []
    for k in 1:mel
        if any(mv[k, 8] .== k4)
            e4 = [e4, k]
        end
    end
    ef4 = 4 * ones(size(e4))

    bound = sort([k1; k2; k3; k4])
    mbound = [e1' ef1'; e2' ef2'; e3' ef3'; e4' ef4']

    ## specify boundary information for graphics
    # bndxy: (x,y)-coordinates of vertices that define the domain and obstacle(s)
    bndxy = [-1 -1; 1 -1; 1 1; -1 1]

    ## centroid coordinate vector
    # xyp = q2p1grid(x, y, xy, mv, bound)

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

    # plotting of the grid
    # adj = spzeros(nvtx, nvtx)
    # for i = 1:nel
    #    adj[convert(Int, mv[i, 1]), convert(Int, mv[i, 2])] = 1
    #    adj[convert(Int, mv[i, 2]), convert(Int, mv[i, 3])] = 1
    #    adj[convert(Int, mv[i, 3]), convert(Int, mv[i, 4])] = 1
    #    adj[convert(Int, mv[i, 4]), convert(Int, mv[i, 1])] = 1
    # end

    # stokes q2-p1 matrix generator
    (ae, bxe, bye, bbxe, bbye) = stokes_q2p1(xy, xyp, mv)

    return (xy, xyp, mv, bound, ae, bxe, bye, bbxe, bbye)

end
