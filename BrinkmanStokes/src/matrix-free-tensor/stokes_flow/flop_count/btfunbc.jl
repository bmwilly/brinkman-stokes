###BTFUNBC
#input
# u           vector to multiply
# xy          Q2 nodal coordinate vector
# xyp         Q1 nodal coordinate vector
# mv          Q2 element mapping matrix
# bound       indices of boundary points
# bxe         local Q2-Q1 divergence matrix
# bye         local Q2-Q1 divergence matrix
#output
# w           B' * u
function btfunbc(u, kparams)

    xy = kparams["xy"]
    xyp = kparams["xyp"]
    mv = kparams["mv"]
    bound = kparams["bound"]
    bxe = kparams["bxe"]
    bye = kparams["bye"]

    # get variables
    nvtx = length(xy[:, 1])
    nu = 2nvtx
    np = 3length(xyp[:, 1])
    nel = length(mv[:, 1])
    mp = [[1:3:3nel] [2:3:3nel] [3:3:3nel]]

    up = u[nu+1:nu+np]
    bxes = squeeze(bxe[1, :, :], 1)
    byes = squeeze(bye[1, :, :], 1)

    w = zeros(nu + np, 1)

    nflops = 0

    # w = SharedArray(Float64, (nu + np, 1), pids = workers())
    # nelworker = nel/nworkers()

    # @parallel for worker = 1:nworkers()

    # for e = (nelworker*worker - nelworker + 1):(nelworker*worker)
    for e = 1:nel
        ind3 = mp[e, :]'
        ind9 = mv[e, :]'
        indbd = findall(in(bound), ind9)

        bxesbd = copy(bxes)
        byesbd = copy(byes)
        nbd = length(indbd)
        bxesbd[:, indbd] = zeros(3, nbd)
        byesbd[:, indbd] = zeros(3, nbd)

        wepx = up[ind3]' * bxesbd
        wepy = up[ind3]' * byesbd

        m, n = size(bxesbd)
        nflops += m * (2n - 1)
        m, n = size(byesbd)
        nflops += m * (2n - 1)

        w[ind9] += wepx'
        w[ind9+nvtx] += wepy'

        nflops += length(wepx')
        nflops += length(wepy')
    end

    # end # parallel

    vec(w), nflops
end
