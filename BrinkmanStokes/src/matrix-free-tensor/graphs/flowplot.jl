using PyPlot

function flowplot(sol, kparams)

  x = kparams["x"]; y = kparams["y"];
  xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"]
  ae = kparams["ae"]; bxe = kparams["bxe"]; bye = kparams["bye"]; bbxe = kparams["bbxe"]; bbye = kparams["bbye"]
  nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])

  u = sol[1:nu]
  p = sol[nu+1:end]
  ux = reshape(u[1:nvtx], length(x), length(y))'
  uy = reshape(u[nvtx+1:end], length(x), length(y))'

  # figure();
  # streamplot(x, y, ux, uy, density = 4, color = ux);
  # axis([-1,1,-1,1])

  figure();
  # pcolor(x, y, kp, cmap = "Greys");
  quiver(x, y, ux, uy, ux, scale = 20);
  axis([-1,1,-1,1]);
end
