###BFUNBC
#input
# u           vector to multiply
# xy          Q2 nodal coordinate vector
# xyp         Q1 nodal coordinate vector
# mv          Q2 element mapping matrix
# bound       indices of boundary points
# bxe         local Q2-Q1 divergence matrix
# bye         local Q2-Q1 divergence matrix
#output
# w           B * u
function bfunbc(u, kparams)

    xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"]
    bxe = kparams["bxe"]; bye = kparams["bye"]

    # get variables
    nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
    nel = length(mv[:, 1])
    bxes = squeeze(bxe[1, :, :], 1)
    byes = squeeze(bye[1, :, :], 1)
    mp = [[1:3:3nel] [2:3:3nel] [3:3:3nel]]
    w = zeros(nu + np)
    # w = SharedArray(Float64, (nvtx, 1), pids = workers())

    # zero dirichlet boundary conditions
    ux = u[1:nvtx]
    uy = u[(nvtx + 1):nu]
    ux[bound] = zeros(length(bound))
    uy[bound] = zeros(length(bound))

    # nelworker = nel/nworkers()
    # @parallel for worker = 1:nworkers()
    # for e = (nelworker*worker - nelworker + 1):(nelworker*worker)
    for e in 1:nel
        ind9 = vec(mv[e, :]')
        ind3 = vec(mp[e, :]') + nu
        w[ind3] += bxes * ux[ind9]
        w[ind3] += byes * uy[ind9]
    end
    # end

    return vec(w)
end
