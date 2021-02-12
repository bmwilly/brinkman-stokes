using Dierckx
using GadFly
reload("stokes_flow/streambc.jl")
reload("grids/findobsXY.jl")

###FLOWPLOT2 plots flow data on general domain
# input
#   sol           flow solution vector
#   By            velocity y-derivative matrix
#   Bx            velocity x-derivative matrix
#   A             vector diffusion matrix
#   xy            velocity nodal coordinate vector
#   xyp           pressure nodal coordinate vector
#   x             vector of x-axis interpolation points
#   y             vector of y-axis interpolation points
#   bound         boundary vertex vector
function flowplot2(sol, domain)

    xst = sol["xst"]; By = sol["By"]; Bx = sol["Bx"]; A = sol["A"];
    xy = sol["xy"]; xyp = sol["xyp"]; x = sol["x"]; y = sol["y"];
    bound = sol["bound"]; bndxy = sol["bndxy"]; bnde = sol["bnde"];
    obs = sol["obs"]

    nvtx = length(xy[:, 1]); nu = 2nvtx;
    Asv = A[1:nvtx, 1:nvtx]

  # compute auxiliary quantities
    u = xst[1:nu]
    p = xst[nu + 1:end]
    f = [By -Bx] * u
    (Asv, fsv) = streambc(Asv, f, xy, bound, domain)
    phi = Asv \ fsv

  # # plot pressure
  # p = p[1:3:end]; xx = x[1:end]; yy = y[1:end]
  #
  # # interpolate to a cartesian product mesh
  # (X, Y) = meshgrid(xx, yy)
  # xysol = griddata(xyp[:, 1], xyp[:, 2], p, X, Y)
  # if size(obs, 1) != 0
  #   (II, JJ) = findobsXY(obs, X, Y, bndxy)
  #   xysol[II, JJ] = NaN
  # end

  ## plot pressure
  # p = p[1:3:end]
  # xpsol = unique(xyp[:, 1]); ypsol = unique(xyp[:, 2])
  # xypsol = reshape(p, length(xpsol), length(ypsol))'
  #
  # plot(x = xpsol, y = ypsol, z = xypsol)

  # plot(x = xy[:, 1], y = xy[:, 2], z = phi, Geom.contour)

  # plot velocity

  # grid = RectangleGrid(x, y)
  # xysol = Float64[]
  # for i = 1:length(xy[:, 1])
  #   pt = interpolate(grid, phi, vec([xy[i, 1], xy[i, 2]]))
  #   push!(xysol, pt)
  # end

  # xsol = unique(x)
  # ysol = unique(y)

    spl = Spline2D(xy[:, 1], xy[:, 2], phi, s=length(xy[:, 1]))
    z = evalgrid(spl, vec(x), vec(y))

  # if size(obs, 1) != 0
  #   KK = findobsXY(obs, X, Y, bndxy)
  #   z[KK] = NaN
  # end

    plot(x=x, y=y, z=z, Geom.contour)

  # (X, Y) = meshgrid(x, y)
  # # xysol = griddata(xy[:, 1], xy[:, 2], phi, X, Y)
  # xsol = unique(xy[:, 1]); ysol = unique(xy[:, 2])
  # phi2 = vcat(zeros(length(xsol)*length(ysol) - length(phi), 1), phi)
  # xysol = reshape(phi2, length(xsol), length(ysol))
  # plot(x = xsol, y = ysol, z = xysol, Geom.contour)

  # spline = Spline2D(xy[:, 1], xy[:, 2], phi)
  # spline = Spline2D(x, y, phi)
  # spline = Spline2D(vec(X), vec(Y), phi)
  # spline = Spline2D(vec(x), vec(y), phi)
  # xysol = evalgrid(spline, X, Y)
  # plot(x = X, y = Y, z = xysol, Geom.Contour)

  # if size(obs, 1) != 0
  #   (II, JJ) = findobsXY(obs, X, Y, bndxy)
  #   xysol(II, JJ) = NaN
  # end
  #
  # plot(x = xsol, y = ysol, z = xysol, Geom.contour)
  # plot(x = xy[:, 1], y = xy[:, 2], z = phi, Geom.contour)

end
