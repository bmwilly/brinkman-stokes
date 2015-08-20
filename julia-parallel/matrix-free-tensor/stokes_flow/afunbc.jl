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
function afunbc(u::SharedArray, kparams)

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
    # t = toc();
    # @show t1;
    # t1 += t;

    # tic()
    @sync begin
      for p in procs()
        @async remotecall_wait(p, loop_elem_chunk!, w, u, aes, mv, nvtx)
      end
      # for p in workers()
      #   @async remotecall_wait(p, loop_elem_chunk!, w, u, Ux, Uy, aes, mv, nvtx)
      # end
    end
    # t = toc();
    # @show t2;
    # t2 += t;

    # tic()
    w[bound] = uu[bound]
    w[bound+nvtx] = uu[bound+nvtx]
    vec(w)
    # t = toc();
    # @show t3;
    # t3 += t;
end

@everywhere function myrange(mv::SharedArray)
  ind = indexpids(mv)
  if ind == 0
    # this worker is not assigned a piece
    return 1:0
  end
  nchunks = length(procs(mv))
  splits = [iround(s) for s in linspace(0, size(mv, 1), nchunks + 1)]
  splits[ind]+1:splits[ind+1]
end

@everywhere function loop_elem!(w::SharedArray, u::SharedArray, aes, mv, nvtx, prange::UnitRange)
  # @show erange
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

@everywhere loop_elem_chunk!(w, u, aes, mv, nvtx) = loop_elem!(w, u, aes, mv, nvtx, myrange(mv))
