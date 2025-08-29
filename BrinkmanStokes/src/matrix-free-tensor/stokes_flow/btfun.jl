function btfun(u, kparams)

    xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"]
    bxe = kparams["bxe"]; bye = kparams["bye"]

    # get variables
    nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
    nel = length(mv[:, 1])
    bxes = share(bxe)
    byes = share(bye)
    mp = [[1:3:3nel] [2:3:3nel] [3:3:3nel]]
    w = SharedArray(Float64, nu + np)

    up = u[(nu + 1):(nu + np)]

    (n, m) = size(mp)
    Up = zeros(n, m)

    for e in 1:nel
        ind3 = mp[e, :]'
        Up[e, :] = up[ind3]'
    end
    Wpx = Up * bxes; Wpy = Up * byes
    for e in 1:nel
        ind9 = mv[e, :]'
        w[ind9] += Wpx[e, :]'
        w[ind9 + nvtx] += Wpy[e, :]'
    end

    return vec(w)
end
