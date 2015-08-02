using PyPlot
# using Gadfly

# reload("stokes_flow/bbfunplot.jl")
# reload("stokes_flow/streambc.jl")
# reload("stokes_flow/afunstreambc.jl")

function flowplot(sol, kparams)

  x = kparams["x"]; y = kparams["y"];
  xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"]
  ae = kparams["ae"]; bxe = kparams["bxe"]; bye = kparams["bye"]; bbxe = kparams["bbxe"]; bbye = kparams["bbye"]
  nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])

  u = sol[1:nu]
  p = sol[nu+1:end]
  ux = reshape(u[1:nvtx], length(x), length(y))'
  uy = reshape(u[nvtx+1:end], length(x), length(y))'

  figure(); PyPlot.streamplot(x, y, ux, uy, density = 2, color = ux); colorbar();

  # x = xy[:, 1]
  # y = xy[:, 2]
  # xp = xyp[:, 1]
  # yp = xyp[:, 2]
  # nvtx = length(x)
  # nu = 2nvtx
  # np = 3length(xp)
  #
  # # compute auxilliary quantities
  # u = sol[1:nu]
  # p = sol[nu + 1:end]
  # f = bbfunplot(u, kparams)
  # fsv = streambc(f, kparams)
  #
  # A = u -> afunstreambc(u, kparams)
  # A = LinearOperator(length(fsv), Float64, A)
  # (phi, flag, err, iter, resvec) = gmres(
  #   A, fsv, length(fsv);
  #   tol = 1e-6, maxIter = 100,
  # )
  #
  # ## plot pressure
  # p = p[1:3:end]
  # xpsol = unique(xp)
  # ypsol = unique(yp)
  # xypsol = reshape(p, length(xpsol), length(ypsol))'
  #
  # # figure()
  # # mesh(xpsol, ypsol, xypsol)
  # # savefig("graphs/pressure.png")
  #
  # ## plot velocity
  # xsol = unique(x)
  # ysol = unique(y)
  # xysol = reshape(phi, length(xsol), length(ysol))'
  #
  # # figure()
  # # contour(xsol, ysol, xysol)
  # # savefig("graphs/velocity.png")
  #
  # plot(z = xysol, Geom.contour)

end
