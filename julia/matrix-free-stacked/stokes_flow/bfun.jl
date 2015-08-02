function bfun(u, kparams)

  xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"];
  bxe = kparams["bxe"]; bye = kparams["bye"]

  # get variables
  nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
  nel = length(mv[:, 1])
  bxes = squeeze(bxe[1, :, :], 1)
  byes = squeeze(bye[1, :, :], 1)
  mp = [[1:3:3nel] [2:3:3nel] [3:3:3nel]]
  w = zeros(nu+np)
  # w = SharedArray(Float64, (nvtx, 1), pids = workers())

  ux = u[1:nvtx]
  uy = u[nvtx + 1:nu]

  # @parallel for e = 1:nel
  for e = 1:nel
      ind9 = vec(mv[e, :]')
      ind3 = vec(mp[e, :]')
      w[ind3] += bxes * ux[ind9]
      w[ind3] += byes * uy[ind9]
  end

  vec(w)
end
