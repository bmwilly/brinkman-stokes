include("../stokes_flow/solve_stokes.jl")
include("../graphs/flowplot.jl")
include("../helpers/input.jl")
using LinearOperators
using KrylovMethods
include("../stokes_flow/square_stokes.jl")
include("../stokes_flow/fbc.jl")
include("../stokes_flow/gbc.jl")
include("../stokes_flow/kfunbc.jl")
include("../solvers/mg_diff.jl")
include("../solvers/m_st_mg.jl")
include("../stokes_flow/qfun.jl")
include("../stokes_flow/qfun_diag.jl")

kparams = square_stokes()
u = linspace(1, 62, 62)

x = kparams["x"]; y = kparams["y"];
xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"]
fst = fbc(kparams)
gst = gbc(kparams)
rhs = vec([fst; gst])
nu = length(xy[:, 1]); nv = 2nu; np = 3length(xyp[:, 1])

# set up GMG preconditioner
(mgdata, smooth_data, sweeps, stype, npre, npost, nc) = mg_diff(kparams)

mparams = {
  "nv" => nv,
  "np" => np,
  "Q" => u -> qfun(u, kparams),
  "mgdata" => mgdata,
  "smooth_data" => smooth_data,
  "nc" => nc,
  "npre" => npre,
  "npost" => npost,
  "sweeps" => sweeps
}
w = m_st_mg(u, mparams)

# w = qfun(u, kparams)


# w = kfunbc(u, kparams)

# wa = afunbc(u, kparams); wa = wa[1:50]
# wb = bfunbc(u, kparams)
# wbt = btfunbc(u, kparams)
# w1 = kfunbc(u, kparams)

# include("../../efficient-operators/stokes_flow/square_stokes.jl")
# include("../../efficient-operators/stokes_flow/obstacle_stokes.jl")
# include("../../efficient-operators/stokes_flow/brinkman_stokes.jl")
# include("../../efficient-operators/stokes_flow/flowbc.jl")
# include("../../efficient-operators/solvers/mg_diff.jl")
# include("../../efficient-operators/solvers/m_st_mg.jl")
#
# mats = square_stokes()
# A = mats["A"]; B = mats["B"]; Bx = mats["Bx"]; By = mats["By"];
#   f = mats["f"]; g = mats["g"]; xy = mats["xy"]; xyp = mats["xyp"];
#   bound = mats["bound"]; x = mats["x"]; y = mats["y"];
#   (Ast, Bst, fst, gst) = flowbc(A, B, f, g, xy, bound, 1)
#   np = length(gst)
#   K = [Ast Bst'; Bst spzeros(np, np)]
#   u = linspace(1, 62, 62)
#   w2 = K*u
#
# find(abs(w1 - w2) .> 1e-5)
