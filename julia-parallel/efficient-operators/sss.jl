# reload("stokes_flow/solve_stokes.jl")
# reload("graphs/flowplot.jl")
# reload("helpers/helper_functions.jl")
include("stokes_flow/solve_stokes.jl")
include("graphs/flowplot.jl")
include("helpers/helper_functions.jl")

domain = int(input("Choose domain (1/lid-driven cavity, 2/brinkman): "))
msize = int(input("Mesh size: "))
sol = solve_stokes(domain, msize)
println("done")
# flowplot(sol, domain)
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

writecsv("$(homedir())/Documents/brinkman-stokes/julia-parallel/efficient-operators/temp/sol/x.csv", x)
writecsv("$(homedir())/Documents/brinkman-stokes/julia-parallel/efficient-operators/temp/sol/y.csv", y)
writecsv("$(homedir())/Documents/brinkman-stokes/julia-parallel/efficient-operators/temp/sol/ux.csv", ux)
writecsv("$(homedir())/Documents/brinkman-stokes/julia-parallel/efficient-operators/temp/sol/uy.csv", uy)
