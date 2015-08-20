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
function afunbc(u, kparams)

    tic()
    xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"];
    ae = kparams["ae"]

    # get variables
    nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
    nel = length(mv[:, 1])
    aes = squeeze(ae[1, :, :], 1)
    w = zeros(nu+np)

    # zero dirichlet boundary conditions
    uu = copy(u)
    u[bound] = zeros(length(bound))
    u[bound+nvtx] = zeros(length(bound))

    n,m = size(mv)
    Ux = zeros(m, n); Uy = zeros(m, n)
    t1 = toc();
    @show t1;

    tic()
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
    t2 = toc();
    @show t2;

    tic()
    w[bound] = uu[bound]
    w[bound+nvtx] = uu[bound+nvtx]
    vec(w)
    t3 = toc();
    @show t3;
end
