function bbfunplot(u, kparams)

  xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"];
  bbxe = kparams["bbxe"]; bbye = kparams["bbye"]

  # get variables
  nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
  nel = length(mv[:, 1])
  bbxes = share(bbxe)
  bbyes = share(bbye)
  w = zeros(nvtx)
  # w = SharedArray(Float64, (nvtx, 1), pids = workers())

  ux = u[1:nvtx]
  uy = u[nvtx + 1:nu]

  n,m = size(mv)
  Ux = zeros(m, n); Uy = zeros(m, n)

  for e = 1:nel
    ind = vec(mv[e, :]')
    Uy[:, e] = ux[ind]
    Ux[:, e] = uy[ind]
  end
  Wx = bbxes * Ux; Wy = bbyes * Uy
  for e = 1:nel
    ind = vec(mv[e, :]')
    w[ind] += Wy[:, e]
    w[ind] -= Wx[:, e]
  end

  vec(w)
end
