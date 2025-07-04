###AFUNBC
# input
# u           vector to multiply
# xy          Q2 nodal coordinate vector
# xyp         Q1 nodal coordinate vector
# mv          Q2 element mapping matrix
# bound       indices of boundary points
# ae          local Q2 diffusion derivative matrix
# output
# w           A * u
function afunbc(u, kparams)

  xy = kparams["xy"]
  xyp = kparams["xyp"]
  mv = kparams["mv"]
  bound = kparams["bound"]
  ae = kparams["ae"]

  # get variables
  nvtx = length(xy[:, 1])
  nu = 2nvtx
  np = 3length(xyp[:, 1])
  nel = length(mv[:, 1])
  aes = squeeze(ae[1, :, :], 1)
  w = zeros(nu + np, 1)
  # w = SharedArray(Float64, (nu + np, 1), pids = workers())
  # nelworker = nel/nworkers()

  nflops = 0

  # @parallel for worker = 1:nworkers()

  #   for e = (nelworker*worker - nelworker + 1):(nelworker*worker)
  # @parallel for e = 1:nel
  # for e = 1:nel/nworkers()
  for e = 1:nel
    ind = mv[e, :]'
    indbd = findall(in(bound), ind)
    indint = setdiff(int(linspace(1, 9, 9)), indbd)
    indb = ind[indbd]

    ux = u[ind]
    ux[indbd] = zeros(1, length(indbd))

    uy = u[ind+nvtx]
    uy[indbd] = zeros(1, length(indbd))

    wex = aes * ux
    wey = aes * uy

    m, n = size(aes)
    nflops += m * (2n - 1)
    m, n = size(aes)
    nflops += m * (2n - 1)

    # impose boundary conditions
    wex[indbd] = u[indb]
    wey[indbd] = u[indb+nvtx]

    # add interior points
    w[ind[indint]] += wex[indint]
    w[ind[indint]+nvtx] += wey[indint]

    nflops += length(wex[indint])
    nflops += length(wey[indint])

    # add boundary points only when entry == 0
    wbx = w[indb]
    wby = w[indb+nvtx]
    ind0x = findall(wbx .== 0)
    ind0y = findall(wby .== 0)
    webx = wex[indbd]
    weby = wey[indbd]
    w[indb[ind0x]] += webx[ind0x]
    w[indb[ind0y]+nvtx] += weby[ind0y]

    nflops += length(webx[ind0x])
    nflops += length(weby[ind0y])
  end


  # end # parallel

  vec(w), nflops
end
