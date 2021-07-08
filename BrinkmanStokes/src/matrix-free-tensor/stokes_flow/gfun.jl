###GFUN
#input
    # u           vector to multiply
    # xy          Q2 nodal coordinate vector
    # xyp         Q1 nodal coordinate vector
    # mv          Q2 element mapping matrix
    # bound       indices of boundary points
    # ge          local Q2 vector mass matrix
#output
    # w           G * u
function gfun(u, kparams)

  xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"]
  ge = kparams["ge"];

  # get variables
  nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
  nel = length(mv[:, 1])
  ges = share(ge)
  w = zeros(nu+np)

  n,m = size(mv)
  Ux = zeros(m, n); Uy = zeros(m, n)

  for e = 1:nel
    ind = vec(mv[e, :]')
    Ux[:, e] = u[ind]
    Uy[:, e] = u[ind+nvtx]
  end
  Wx = ges * Ux; Wy = ges * Uy
  for e = 1:nel
    ind = vec(mv[e, :]')
    w[ind] += Wx[:, e]
    w[ind+nvtx] += Wy[:, e]
  end

  vec(w)
end
