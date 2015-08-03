function btfun(u, kparams)

    xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"];
    bxe = kparams["bxe"]; bye = kparams["bye"]

    # get variables
    nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
    nel = length(mv[:, 1])
    bxes = squeeze(bxe[1, :, :], 1)
    byes = squeeze(bye[1, :, :], 1)
    mp = [[1:3:3nel] [2:3:3nel] [3:3:3nel]]
    w = zeros(nu+np)
    # w = SharedArray(Float64, (nu + np, 1), pids = workers())

    up = u[nu + 1:nu + np]

    # @parallel for e = 1:nel
    for e = 1:nel
        ind3 = vec(mp[e, :]')
        ind9 = vec(mv[e, :]')

        wepx = up[ind3]' * bxes;
        wepy = up[ind3]' * byes;
        w[ind9] += wepx'
        w[ind9 + nvtx] += wepy'
    end

    vec(w)
end
