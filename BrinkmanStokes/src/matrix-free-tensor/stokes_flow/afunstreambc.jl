###AFUNSTREAMBC
#input
# u           vector to multiply
# xy          Q2 nodal coordinate vector
# xyp         Q1 nodal coordinate vector
# mv          Q2 element mapping matrix
# bound       indices of boundary points
# ae          local Q2 diffusion derivative matrix
#output
# w           A * u
function afunstreambc(u, kparams)

    xy = kparams["xy"]
    xyp = kparams["xyp"]
    mv = share(kparams["mv"])
    bound = kparams["bound"]
    ae = kparams["ae"]

    # get variables
    nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
    nel = length(mv[:, 1])
    aes = share(ae)
    w = SharedArray(Float64, nu + np)

    # zero dirichlet boundary conditions
    uu = copy(u)
    u[bound] = zeros(length(bound))

    @sync begin
        for p in procs()
            @async remotecall_wait(p, afunstreambc_loop_chunk!, w, u, aes, mv)
        end
    end

    w[bound] = uu[bound]
    return vec(w)
end

@everywhere function afunstreambc_loop!(w, u, aes, mv, prange)
    pnel = length(prange)
    mvp = mv[prange, :]
    U = zeros(size(mvp, 2), pnel)
    for e in 1:pnel
        ind = vec(mvp[e, :]')
        U[:, e] = u[ind]
    end
    W = aes * U
    for e in 1:pnel
        ind = vec(mvp[e, :]')
        w[ind] += W[:, e]
    end
    w
end

@everywhere afunstreambc_loop_chunk!(w, u, aes, mv) = afunstreambc_loop!(w, u, aes, mv, myrange(mv))
