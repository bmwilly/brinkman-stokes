using PyPlot
# using Dierckx
# using Gaston
# using Gadfly
# using Plotly

reload("stokes_flow/streambc.jl")

function flowplot(sol, domain)

  xst = sol["xst"]; By = sol["By"]; Bx = sol["Bx"]; A = sol["A"];
  xy = sol["xy"]; xyp = sol["xyp"]; x = sol["x"]; y = sol["y"];
  bound = sol["bound"]; bndxy = sol["bndxy"]; bnde = sol["bnde"];
  obs = sol["obs"];

  nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1]);
  Asv = A[1:nvtx, 1:nvtx]; x = vec(x); y = vec(y);
  xp = unique(xyp[:, 1]); yp = unique(xyp[:, 2]);

  # compute auxilliary quantities
  u = xst[1:nu];
  p = xst[nu + 1:end];
  # f = [By -Bx] * u
  # (Asv, fsv) = streambc(Asv, f, xy, bound, domain)
  # phi = Asv\fsv

  ## plot pressure
  # p = p[1:3:end]
  # zp = reshape(p, length(xp), length(yp))'
  # p1 = Gadfly.plot(x = xp, y = yp, z = zp, Geom.contour);

  ## plot velocity
  ux = reshape(u[1:nvtx], length(x), length(y))';
  uy = reshape(u[nvtx+1:end], length(x), length(y))';
  figure(); streamplot(x, y, ux, uy, density = 2, color = ux); colorbar()

  ## plot velocity
  # ux = reshape(u[1:nvtx], length(x), length(y))
  # p2 = plot(x = x, y = y, z = ux, Geom.contour);

  # X,Y = meshgrid(x, y)
  # PyPlot.surf(X, Y, ux)
  # PyPlot.contourf(X,Y,ux)

  # data = [[
  #   "x" => x,
  #   "y" => y,
  #   "z" => ux,
  #   "type" => "contour"
  # ]]
  # r1 = Plotly.plot(data)
  # url1 = r1["url"]
  #
  # data = [[
  #   "z" => ux,
  #   "type" => "surface"
  # ]]
  # r2 = Plotly.plot(data)
  # url2 = r2["url"]
  #
  # uy = reshape(u[nvtx+1:end], length(x), length(y))'
  # p3 = plot(x = x, y = y, z = uy, Geom.contour);
  #
  # z = reshape(phi, length(x), length(y))'
  # p4 = plot(x = x, y = y, z = z, Geom.contour);

  if domain == 1
    ## plot pressure
    # p = p[1:3:end]
    # zp = reshape(p, length(xp), length(yp))'
    # p1 = plot(x = xp, y = yp, z = zp, Geom.contour);
    # draw(PNG("graphs/cavity_pressure.png", 8inch, 4inch), p1);

    ## plot velocity
    # ux = reshape(u[1:nvtx], length(x), length(y))'
    # p2 = plot(x = x, y = y, z = ux, Geom.contour);
    # draw(PNG("graphs/cavity_velocity_x.png", 8inch, 4inch), p2);

    # uy = reshape(u[nvtx+1:end], length(x), length(y))'
    # p3 = plot(x = x, y = y, z = uy, Geom.contour);
    # draw(PNG("graphs/cavity_velocity_y.png", 8inch, 4inch), p3);

    # ## plot pressure
    # p = p[1:3:end]
    # zp = reshape(p, length(xp), length(yp))'
    # p1 = plot(x = xp, y = yp, z = zp, Geom.contour)
    # draw(PNG("graphs/cavity_pressure.png", 8inch, 4inch), p1)
    #
    # ## plot velocity
    # z = reshape(phi, length(x), length(y))'
    # p2 = plot(x = x, y = y, z = z, Geom.contour)
    # # p2 = plot(x = x, y = y, z = z, Stat.contour(levels = 50))
    # draw(PNG("graphs/cavity_velocity.png", 8inch, 4inch), p2)
  elseif domain == 2
    # plot velocity
    # spl = Spline2D(
    #   xy[:, 1], xy[:, 2], phi;
    #   kx = 4, ky = 4, s = length(xy[:, 1])
    # )
    # z = evalgrid(spl, x, y)
    #
    # (X, Y) = meshgrid(x, y)
    # if size(obs, 1) != 0
    #   KK = findobsXY(obs, X, Y, bndxy)
    #   phi[KK] = NaN
    # end
    # p2 = plot(
    #   layer(x = x, y = y, z = z, Geom.contour),
    #   layer(x = [bndxy[5:8, 1], bndxy[5, 1]], y = [bndxy[5:8, 2], bndxy[5, 2]], Geom.path)
    # )
    # p2 = plot(x = x, y = y, z = z, Geom.contour)
    # draw(PNG("graphs/obstacle_velocity.png", 8inch, 4inch), p2)
  elseif domain == 3
    ## plot pressure
    # p = p[1:3:end]
    # zp = reshape(p, length(xp), length(yp))'
    # p1 = plot(x = xp, y = yp, z = zp, Geom.contour);
    # draw(PNG("graphs/brinkman_pressure.png", 8inch, 4inch), p1);
    # draw(PNG("graphs/brinkman_streams.png", 8inch, 4inch), p4);

    P = sol["P"];
    # figure(); pcolor(P')

    # data = [["x" => x, "y" => y, "z" => ux, "type" => "contour"]]
    # rx = Plotly.plot(data)
    # urlx = rx["url"]
    #
    # data = [[
    #   "x" => x,
    #   "y" => y,
    #   "z" => z,
    #   "type" => "contour"
    # ]]
    # r1 = Plotly.plot(data)
    # url1 = r1["url"]
    #
    # data = [[
    #   "z" => z,
    #   "type" => "surface"
    # ]]
    # r2 = Plotly.plot(data)
    # url2 = r2["url"]
    #
    # return(urlx)
  end
end
