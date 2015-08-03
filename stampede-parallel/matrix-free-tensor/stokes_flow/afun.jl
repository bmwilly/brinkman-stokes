function afun(u, kparams)

  xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"];
  ae = kparams["ae"]

  # get variables
  nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
  nel = length(mv[:, 1])
  aes = squeeze(ae[1, :, :], 1)
  w = zeros(nu+np)
  # w = SharedArray(Float64, (nu + np, 1), pids = workers())

  n,m = size(mv)
  Ux = zeros(m, n); Uy = zeros(m, n)

  # @parallel for e = 1:nel
  # @spawn for e = 1:nel
  for e = 1:nel
    ind = vec(mv[e, :]')
    Ux[:, e] = u[ind]
    Uy[:, e] = u[ind+nvtx]
  end

  Wx = aes * Ux; Wy = aes * Uy

  for e = 1:nel
    ind = vec(mv[e, :]')
    w[ind] += Wx[:, e]
    w[ind+nvtx] += Wy[:, e]
  end

  vec(w)
end
