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

    xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"];
    ae = kparams["ae"]

    # get variables
    nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
    nel = length(mv[:, 1])
    aes = squeeze(ae[1, :, :], 1)
    w = zeros(nu+np)

    # zero dirichlet boundary conditions
    uu = copy(u)
    u[bound] = zeros(length(bound))
    u[bound+nvtx] = zeros(length(bound))
    # w = SharedArray(Float64, (nu + np, 1), pids = workers())
    # nelworker = nel/nworkers()

    # @parallel for worker = 1:nworkers()

    #   for e = (nelworker*worker - nelworker + 1):(nelworker*worker)
      # @parallel for e = 1:nel
      # for e = 1:nel/nworkers()
    for e = 1:nel
      ind = vec(mv[e, :]')
      w[ind] += aes * u[ind]
      w[ind+nvtx] += aes * u[ind+nvtx]
    end

    w[bound] = uu[bound]
    w[bound+nvtx] = uu[bound+nvtx]
    # end # parallel

    vec(w)
end
