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
  ges = squeeze(ge[1, :, :], 1)
  w = zeros(nu+np)

  for e = 1:nel
    ind = vec(mv[e, :]')
    w[ind] += ges * u[ind]
    w[ind+nvtx] += ges * u[ind+nvtx]
  end

  vec(w)
end
