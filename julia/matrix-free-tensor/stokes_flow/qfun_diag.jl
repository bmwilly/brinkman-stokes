###QFUN_DIAG
#input
    # u           vector to multiply
    # xy          Q2 nodal coordinate vector
    # xyp         Q1 nodal coordinate vector
    # mv          Q2 element mapping matrix
    # bound       indices of boundary points
    # qe          local Q1 mass matrix
#output
    # w           Q * u
function qfun_diag(u, kparams)

  xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"]
  qe = kparams["qe"]

  nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
  nel = length(mv[:, 1])
  mp = [[1:3:3nel] [2:3:3nel] [3:3:3nel]]
  qes = squeeze(qe[1, :, :], 1)
  w = zeros(nu+np)

  n,m = size(mv)
  U = zeros(m, n)

  for e = 1:nel
    ind = vec(mp[e, :]') + nu
    U[:, e] = u[ind]
  end
  W = diagm(diag(qes)) * U
  for e = 1:nel
    ind = vec(mp[e, :]') + nu
    w[ind] += W[:, e]
  end

  w = vec(w[nu+1:nu+np])
end
