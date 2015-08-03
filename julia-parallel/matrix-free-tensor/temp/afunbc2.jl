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
function afunbc2(u, kparams)

    xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"];
    ae = kparams["ae"]

    # get variables
    nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
    nel = length(mv[:, 1])
    aes = squeeze(ae[1, :, :], 1)
    # w = zeros(nu+np)
    w = zeros(nvtx)

    # zero dirichlet boundary conditions
    u = u[1:nvtx]
    uu = copy(u)
    u[bound] = zeros(length(bound))
    # u[bound+nvtx] = zeros(length(bound))

    n,m = size(mv)
    U = zeros(m, n)

    for e = 1:nel
      ind = vec(mv[e, :]')
      U[:, e] = u[ind]
    end

    W = aes * U

    for e = 1:nel
      ind = vec(mv[e, :]')
      w[ind] += W[:, e]
    end

    w[bound] = uu[bound]

    vec(w)
end

function afunbc3(u, kparams)
  xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"];
  ae = kparams["ae"]

  # get variables
  nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
  nel = length(mv[:, 1])
  aes = squeeze(ae[1, :, :], 1)

  w = afunbc2(u, kparams)
  wout = [w; w; zeros(np)]
  vec(wout)
end
