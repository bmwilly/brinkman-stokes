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

    xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"];
    bxe = kparams["bxe"]; bye = kparams["bye"]

    # get variables
    nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
    nel = length(mv[:, 1])
    mp = [[1:3:3nel] [2:3:3nel] [3:3:3nel]]

    ux = u[1:nu/2]
    uy = u[nu/2 + 1:nu]
    bxes = squeeze(bxe[1, :, :], 1)
    byes = squeeze(bye[1, :, :], 1)

    w = zeros(nu + np, 1)
    # w = SharedArray(Float64, (nu + np, 1), pids = workers())
    # nelworker = nel/nworkers()

    nflops = 0

    # @parallel for worker = 1:nworkers()


    # for e = (nelworker*worker - nelworker + 1):(nelworker*worker)
    for e = 1:nel
        ind9 = mv[e, :]'
        indbd = findin(ind9, bound)
        indint = setdiff(int(linspace(1, 9, 9)), indbd)

        bxesbd = copy(bxes)
        byesbd = copy(byes)
        nbd = length(indbd)
        bxesbd[:, indbd] = zeros(3, nbd)
        byesbd[:, indbd] = zeros(3, nbd)
        wex = bxesbd * ux[ind9]
        wey = byesbd * uy[ind9]

        m,n = size(bxesbd); nflops += m*(2n-1)
        m,n = size(byesbd); nflops += m*(2n-1)


        ind3 = mp[e, :]' + nu
        w[ind3] += wex
        w[ind3] += wey

        nflops += length(wex)
        nflops += length(wey)
    end

    # end

    vec(w), nflops
end
