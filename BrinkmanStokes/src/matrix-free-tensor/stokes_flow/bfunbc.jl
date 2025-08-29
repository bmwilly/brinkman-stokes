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
    bxes = share(bxe)
    byes = share(bye)
    mp = [[1:3:3nel] [2:3:3nel] [3:3:3nel]]
    w = zeros(nu + np)
    # w = SharedArray(Float64, (nvtx, 1), pids = workers())

    # zero dirichlet boundary conditions
    ux = u[1:nvtx]
    uy = u[(nvtx + 1):nu]
    ux[bound] = zeros(length(bound))
    uy[bound] = zeros(length(bound))

    n, m = size(mv)
    Ux = zeros(m, n); Uy = zeros(m, n)

    for e in 1:nel
        ind9 = vec(mv[e, :]')
        Ux[:, e] = ux[ind9]
        Uy[:, e] = uy[ind9]
    end
    Wx = bxes * Ux; Wy = byes * Uy
    for e in 1:nel
        ind3 = vec(mp[e, :]') + nu
        w[ind3] += Wx[:, e]
        w[ind3] += Wy[:, e]
    end

    return vec(w)
end
