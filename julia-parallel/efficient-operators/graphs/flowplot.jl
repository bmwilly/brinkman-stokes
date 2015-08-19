using PyPlot
# using Dierckx
# using Gaston
# using Gadfly
# using Plotly

function flowplot(sol, domain)

  xst = sol["xst"]; By = sol["By"]; Bx = sol["Bx"]; A = sol["A"];
  xy = sol["xy"]; xyp = sol["xyp"]; x = sol["x"]; y = sol["y"];
  bound = sol["bound"]; bndxy = sol["bndxy"]; bnde = sol["bnde"];
  obs = sol["obs"];
  if domain == 3; kappa = sol["kappa"]; kp = reshape(kappa, length(x), length(y))'; end;

  nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1]);
  Asv = A[1:nvtx, 1:nvtx]; x = vec(x); y = vec(y);
  xp = unique(xyp[:, 1]); yp = unique(xyp[:, 2]);

  # compute auxilliary quantities
  u = xst[1:nu];
  p = xst[nu + 1:end];

  ## plot pressure
  # p = p[1:3:end]
  # zp = reshape(p, length(xp), length(yp))'
  # figure(); streamplot(xp, yp, zp)

  ## plot velocity
  ux = reshape(u[1:nvtx], length(x), length(y))';
  uy = reshape(u[nvtx+1:end], length(x), length(y))';

  # writecsv("home/bmw313/Documents/tmp/sol/x.csv", x)
  # writecsv("home/bmw313/Documents/tmp/sol/y.csv", y)
  # writecsv("home/bmw313/Documents/tmp/sol/ux.csv", ux)
  # writecsv("home/bmw313/Documents/tmp/sol/uy.csv", uy)


  figure(); streamplot(x, y, ux, uy, density = 4, color = ux); axis([-1,1,-1,1]);
  # figure(); streamplot(x, y, ux, uy, density = 4, color = ux); quiver(x, y, ux, uy); axis([-1,1,-1,1])
  # figure(); quiver(x, y, ux, uy, scale = 20); axis([-1,1,-1,1])
  # figure(); streamplot(x, y, ux, uy, density = 4); quiver(x, y, ux, uy, ux); axis([-1,1,-1,1])

  figure();
  # pcolor(x, y, kp, cmap = "Greys");
  quiver(x, y, ux, uy, ux, scale = 2);
  axis([-1,1,-1,1]);
end
