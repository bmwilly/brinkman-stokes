reload("stokes_flow/bfun.jl")

function gbc(kparams)

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
    bcy = 0 * ybd
    k = find((ybd .== 1) & (xbd .> -1) & (xbd .< 1))
    bcx[k] = (1 - xbd[k] .* xbd[k]) .* (1 + xbd[k] .* xbd[k])

    # impose boundary conditions
    bccx = zeros(nvtx)
    bccx[bound] = bcx
    bccy = zeros(nvtx)
    bccy[bound] = bcy

    bc = [bccx; bccy; zeros(np, 1)]
    w -= bfun(bc, kparams)
    w = vec(w[nu+1:nu+np])
end
