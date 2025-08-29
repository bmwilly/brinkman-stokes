###AFUNBC
#input
# u           vector to multiply
# xy          Q2 nodal coordinate vector
# xyp         Q1 nodal coordinate vector
# mv          Q2 element mapping matrix
# bound       indices of boundary points
# ae          local Q2 diffusion derivative matrix
#output
# w           A * u
function afunbc(u, kparams)

    xy = kparams["xy"]; xyp = kparams["xyp"]
    mv = share(kparams["mv"])
    bound = kparams["bound"]
    ae = kparams["ae"]

    # get variables
    nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
    nel = length(mv[:, 1])
    aes = squeeze(ae[1, :, :], 1)
    w = SharedArray(Float64, length(u))

    # zero dirichlet boundary conditions
    uu = copy(u)
    u[bound] = zeros(length(bound))
    u[bound + nvtx] = zeros(length(bound))

    # @time (
    @sync begin
        for p in procs()
            @async remotecall_wait(p, loop_elem_chunk!, w, u, aes, mv, nvtx)
        end
    end
    # )

    w[bound] = uu[bound]
    w[bound + nvtx] = uu[bound + nvtx]
    return vec(w)
end

@everywhere function myrange(mv::SharedArray)
    ind = indexpids(mv)
    if ind == 0
        # this worker is not assigned a piece
        return 1:0
    end
    nchunks = length(procs(mv))
    splits = [iround(s) for s in linspace(0, size(mv, 1), nchunks + 1)]
    (splits[ind] + 1):splits[ind + 1]
end

@everywhere function loop_elem!(w, u, aes, mv, nvtx, erange)
    # @show erange
    for e in erange
        ind = vec(mv[e, :]')
        w[ind] += aes * u[ind]
        w[ind + nvtx] += aes * u[ind + nvtx]
    end
    w
end

@everywhere loop_elem_chunk!(w, u, aes, mv, nvtx) = loop_elem!(w, u, aes, mv, nvtx, myrange(mv))
