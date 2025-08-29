function afun(u, kparams)

    xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"]
    ae = kparams["ae"]

    # get variables
    nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
    nel = length(mv[:, 1])
    aes = squeeze(ae[1, :, :], 1)
    w = zeros(nu + np)
    # w = SharedArray(Float64, (nu + np, 1), pids = workers())

    # @parallel for e = 1:nel
    # @spawn for e = 1:nel
    for e in 1:nel
        ind = vec(mv[e, :]')
        we = aes * u[ind]
        w[ind] += we
        w[ind + nvtx] += we
    end

    return vec(w)
end
