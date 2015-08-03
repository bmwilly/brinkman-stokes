reload("stokes_flow/afun.jl")
reload("stokes_flow/btfun.jl")

function fbc(kparams)

    xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"]

    # get variables
    nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
    nel = length(mv[:, 1])
    w = zeros(nu + np, 1)

    # get boundary
    xbd = xy[bound, 1]
    ybd = xy[bound, 2]

    # a regularized cavity
    bcx = 0 * xbd
    bcy = 0 * xbd
    k = find((ybd .== 1) & (xbd .> -1) & (xbd .< 1))
    bcx[k] = (1 - xbd[k] .* xbd[k]) .* (1 + xbd[k] .* xbd[k])

    # impose boundary conditions
    bccx = zeros(nvtx)
    bccx[bound] = bcx
    bccy = zeros(nvtx)
    bccy[bound] = bcy

    bc = [bccx; bccy; zeros(np, 1)]
    w -= afun(bc, kparams)
    w -= btfun(bc, kparams)

    wx = w[1:nvtx]
    wy = w[nvtx + 1:nu]
    wp = w[nu + 1:nu + np]

    wx[bound] = bcx
    wy[bound] = bcy

    w = vec([wx; wy])
end
