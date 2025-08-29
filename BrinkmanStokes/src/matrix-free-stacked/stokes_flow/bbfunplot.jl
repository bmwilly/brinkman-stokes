function bbfunplot(u, kparams)

    xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"]
    bbxe = kparams["bbxe"]; bbye = kparams["bbye"]

    # get variables
    nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
    nel = length(mv[:, 1])
    bbxes = squeeze(bbxe[1, :, :], 1)
    bbyes = squeeze(bbye[1, :, :], 1)
    w = zeros(nvtx)
    # w = SharedArray(Float64, (nvtx, 1), pids = workers())

    ux = u[1:nvtx]
    uy = u[(nvtx + 1):nu]

    # @parallel for e = 1:nel
    for e in 1:nel
        ind = vec(mv[e, :]')
        w[ind] += bbxes * uy[ind]
        w[ind] -= bbyes * ux[ind]
    end

    return vec(w)
end
