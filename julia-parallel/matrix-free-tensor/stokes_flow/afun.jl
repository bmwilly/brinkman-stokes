function afun(u, kparams)

  xy = kparams["xy"];
  xyp = kparams["xyp"];
  mv = share(kparams["mv"]);
  bound = kparams["bound"];
  ae = kparams["ae"]

  # get variables
  nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
  nel = length(mv[:, 1])
  aes = share(ae)
  w = SharedArray(Float64, nu+np)

  @sync begin
    for p in procs()
      @async remotecall_wait(p, afun_loop_chunk!, w, u, aes, mv, nvtx)
    end
  end

  vec(w)
end

@everywhere function afun_loop!(w::SharedArray, u::SharedArray, aes, mv::SharedArray, nvtx, prange::UnitRange)
  pnel = length(prange)
  mve = mv[prange,:]
  Ux = zeros(size(mve,2), pnel)
  Uy = zeros(size(mve,2), pnel)
  for e = 1:pnel
    ind = vec(mve[e,:]')
    Ux[:,e] = u[ind]
    Uy[:,e] = u[ind+nvtx]
  end
  Wx = aes * Ux
  Wy = aes * Uy
  for e = 1:pnel
    ind = vec(mve[e,:]')
    w[ind] += Wx[:,e]
    w[ind+nvtx] += Wy[:,e]
  end
  w
end

@everywhere afun_loop_chunk!(w, u, aes, mv, nvtx) = afun_loop!(w, u, aes, mv, nvtx, myrange(mv))
