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

    # tic()
    xy = kparams["xy"];
    xyp = kparams["xyp"];
    mv = share(kparams["mv"]);
    bound = kparams["bound"];
    ae = kparams["ae"]

    # get variables
    nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
    nel = length(mv[:, 1])
    aes = share(ae);

    # zero dirichlet boundary conditions
    uu = copy(u)
    u[bound] = zeros(length(bound))
    u[bound+nvtx] = zeros(length(bound))

    w = SharedArray(Float64, nu+np)

    @sync begin
      for p in procs()
        @async remotecall_wait(p, afunbc_loop_chunk!, w, u, aes, mv, nvtx)
      end
    end

    w[bound] = uu[bound]
    w[bound+nvtx] = uu[bound+nvtx]
    vec(w)
end

@everywhere function afunbc_loop!(w::SharedArray, u::SharedArray, aes, mv, nvtx, prange::UnitRange)
  pnel = length(prange)
  mve = mv[prange,:]
  Ux = zeros(size(mve,2), pnel)
  Uy = zeros(size(mve,2), pnel)
  for e = 1:pnel
    ind = vec(mve[e, :]')
    Ux[:, e] = u[ind]
    Uy[:, e] = u[ind+nvtx]
  end
  Wx = aes * Ux
  Wy = aes * Uy
  for e = 1:pnel
    ind = vec(mve[e, :]')
    w[ind] += Wx[:, e]
    w[ind+nvtx] += Wy[:, e]
  end
  w
end

@everywhere afunbc_loop_chunk!(w, u, aes, mv, nvtx) = afunbc_loop!(w, u, aes, mv, nvtx, myrange(mv))
