using PyPlot

function flowplot(sol, domain)

    xst = sol["xst"]
    By = sol["By"]
    Bx = sol["Bx"]
    A = sol["A"]
    xy = sol["xy"]
    xyp = sol["xyp"]
    x = Array(sol["x"])
    y = Array(sol["y"])
    bound = sol["bound"]
    bndxy = sol["bndxy"]
    bnde = sol["bnde"]
    obs = sol["obs"]
    if domain == 3
        kappa = sol["kappa"]
        kp = reshape(kappa, length(x), length(y))'
    end

    nvtx = length(xy[:, 1])
    nu = 2nvtx
    np = 3length(xyp[:, 1])
    Asv = A[1:nvtx, 1:nvtx]
    x = vec(x)
    y = vec(y)
    xp = unique(xyp[:, 1])
    yp = unique(xyp[:, 2])

    # compute auxilliary quantities
    u = xst[1:nu]
    p = xst[nu+1:end]

    ## plot velocity
    ux = reshape(u[1:nvtx], length(x), length(y))'
    uy = reshape(u[nvtx+1:end], length(x), length(y))'

    figure()
    streamplot(x, y, ux, uy, density = 4, color = ux)
    axis([-1, 1, -1, 1])

    figure()
    # pcolor(x, y, kp, cmap = "Greys");
    quiver(x, y, ux, uy, ux, scale = 20)
    axis([-1, 1, -1, 1])
end
