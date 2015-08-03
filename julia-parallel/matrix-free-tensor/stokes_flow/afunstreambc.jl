###AFUNSTREAMBC
#input
    # u           vector to multiply
    # xy          Q2 nodal coordinate vector
    # xyp         Q1 nodal coordinate vector
    # mv          Q2 element mapping matrix
    # bound       indices of boundary points
    # ae          local Q2 diffusion derivative matrix
#output
    # w           A * u
function afunstreambc(u, kparams)

  xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"];
  ae = kparams["ae"]

  # get variables
  nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
  nel = length(mv[:, 1])
  aes = squeeze(ae[1, :, :], 1)
  w = zeros(nu+np)
  # w = SharedArray(Float64, (nvtx, 1), pids = workers())

  # zero dirichlet boundary conditions
  uu = copy(u)
  u[bound] = zeros(length(bound))

  n,m = size(mv)
  U = zeros(m, n)

  # @parallel for e = 1:nel
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
