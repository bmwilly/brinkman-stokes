using MAT
using Gadfly
reload("stokes_flow/solve_stokes.jl")
reload("graphs/flowplot.jl")
reload("graphs/flowplot2.jl")
reload("helpers/input.jl")
reload("helpers/meshgrid.jl")
reload("stokes_flow/stokes_brinkman_q2p1.jl")
reload("helpers/input.jl")
reload("grids/q2p1grid.jl")
reload("grids/cavity_domain.jl")
reload("grids/channel_domain.jl")

## sol = solve_stokes(domain)
domain = 3;

## mats = brinkman_stokes()

# generate Q2 grid for square channel
channel_grid = channel_domain()
grid = q2p1grid(channel_grid)

x = grid["x"]; y = grid["y"];
spe10 = matread("data/spe10.mat")
KU = spe10["KU"]; pU = spe10["pU"];
layer = 1;
K = KU[1, 1:((length(x) - 1) / 2), 1:((length(y) - 1) / 2), layer];
K = squeeze(K[1,:,:], 1);
K = K.^(-1);
L = zeros(size(K));
l1, l2 = size(L)
L[l1 / 2 - 1:l1 / 2 + 2, l2 / 2 - 1:l2 / 2 + 2] = repeat([1e4], outer=[4, 4])

# stokes q2-p1 matrix generator
stokes_grid = merge(grid, {"mv" => channel_grid["mv"]})
stokes_mats = stokes_brinkman_q2p1(stokes_grid, L)

bounds = {
  "bound" => channel_grid["bound"],
  "bndxy" => channel_grid["bndxy"],
  "bnde" => channel_grid["bnde"],
  "obs" => channel_grid["obs"]
}

# keys(mats) =
# {"A", "B", "G", "Q", "Bx", "By", "f", "g", "x", "y", "xyp", "bound"}
mats = merge(stokes_mats, grid, bounds)

A = mats["A"]; B = mats["B"]; Bx = mats["Bx"]; By = mats["By"];
f = mats["f"]; g = mats["g"]; xy = mats["xy"]; xyp = mats["xyp"];
bound = mats["bound"]; x = mats["x"]; y = mats["y"];
Q = mats["Q"]

# boundary conditions
println("imposing (enclosed flow) boundary conditions ...")
(Ast, Bst, fst, gst) = flowbc(A, B, f, g, xy, bound, domain)
np = length(gst)
rhs = vec([fst; gst])

## compute solution
K = [Ast Bst'; Bst spzeros(np, np)]

nv = size(Ast, 1); np = size(Q, 1); nu = nv / 2
Agal = Ast[1:nu, 1:nu]
(mgdata, smooth_data, sweeps, stype, npre, npost, nc) = mg_diff(x, y, Agal)

mparams = {
  "nv" => nv,
  "Q" => Q,
  "mgdata" => mgdata,
  "smooth_data" => smooth_data,
  "nc" => nc,
  "npre" => npre,
  "npost" => npost,
  "sweeps" => sweeps
}
# block GMG preconditioner
M = u -> m_st_mg(u, mparams)

tol = 1e-6; maxit = 100;
(xst, flag, err, iter, resvec) = gmres(
  K, rhs, length(rhs);
  tol=tol, maxIter=maxit, M=M, out=1
)

if flag == 0
    println("GMRES reached desired tolerance at iteration $(length(resvec))")
end

p = plot(
  x=1:1:length(resvec), y=log(10, resvec), Geom.line,
  Guide.xlabel("Iteration"), Guide.ylabel("log_10(residual)")
)
draw(PNG("graphs/iters.png", 12inch, 6inch), p)

sol = {
  "K" => K, "Ast" => Ast, "Bst" => Bst, "M" => M,
  "rhs" => rhs, "fst" => fst, "gst" => gst, "xst" => xst
}

sol = merge(mats, sol)

xst = sol["xst"]; By = sol["By"]; Bx = sol["Bx"]; A = sol["A"];
xy = sol["xy"]; xyp = sol["xyp"]; x = sol["x"]; y = sol["y"];
bound = sol["bound"]; bndxy = sol["bndxy"]; bnde = sol["bnde"];
obs = sol["obs"]

nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
Asv = A[1:nvtx, 1:nvtx]; x = vec(x); y = vec(y);
xp = unique(xyp[:, 1]); yp = unique(xyp[:, 2]);

# compute auxilliary quantities
u = xst[1:nu]
p = xst[nu + 1:end]

## plot pressure
p = p[1:3:end]
zp = reshape(p, length(xp), length(yp))'
p1 = plot(x=xp, y=yp, z=zp, Geom.contour);
draw(PNG("graphs/brinkman_pressure.png", 8inch, 4inch), p1);

## plot velocity
ux = reshape(u[1:nvtx], length(x), length(y))'
p2 = plot(x=x, y=y, z=ux, Geom.contour);
draw(PNG("graphs/brinkman_velocity_x.png", 8inch, 4inch), p2);

uy = reshape(u[nvtx + 1:end], length(x), length(y))'
p3 = plot(x=x, y=y, z=uy, Geom.contour);
draw(PNG("graphs/brinkman_velocity_y.png", 8inch, 4inch), p3);

DelimitedFiles.writedlm("/Users/bwilliams/Documents/temp/ux_julia.csv", ux)
DelimitedFiles.writedlm("/Users/bwilliams/Documents/temp/uy_julia.csv", uy)
