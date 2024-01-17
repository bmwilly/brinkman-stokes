using LinearOperators
using KrylovMethods

reload("stokes_flow/square_stokes.jl")
reload("stokes_flow/fbc.jl")
reload("stokes_flow/gbc.jl")
reload("stokes_flow/kfunbc.jl")
reload("stokes_flow/agal.jl")
reload("solvers/mg_diff.jl")
reload("stokes_flow/qfun.jl")
reload("stokes_flow/qfun_diag.jl")
reload("solvers/m_st_mg.jl")

###SOLVE_STOKES solve stokes problem
function solve_stokes_mg()

  kparams = square_stokes()

  x = kparams["x"]; y = kparams["y"];
  xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"]

  println("imposing (enclosed flow) boundary conditions ...")

  # get RHS vector
  fst = fbc(kparams)
  gst = gbc(kparams)
  rhs = vec([fst; gst])

  ## STOKES problem

  # set parameters
  # hard-coded for now
  tol = 1e-6
  maxit = 100

  # specify matrix-free galerkin
  K = u -> kfunbc(u, kparams)

  # construct preconditioner
  println("multigrid preconditioning...")
  nu = length(xy[:, 1]); nv = 2nu; np = 3length(xyp[:, 1])

  # set up GMG preconditioner
  (mgdata, smooth_data, sweeps, stype, npre, npost, nc) = mg_diff(kparams)

  mparams = {
    "nv" => nv,
    "np" => np,
    "Q" => u -> qfun(u, kparams),
    "Qdiag" => u -> qfun_diag(u, kparams),
    "mgdata" => mgdata,
    "smooth_data" => smooth_data,
    "nc" => nc,
    "npre" => npre,
    "npost" => npost,
    "sweeps" => sweeps
  }

  println("GMRES iteration...")

  # zero initial guess
  x0 = vec(zeros(size(rhs)))

  # MG preconditioner
  M = u -> m_st_mg(u, mparams)

  K = LinearOperator(nv+np, Float64, K)
  M = LinearOperator(nv+np, Float64, M)

  (xst, flag, err, iter, resvec) = gmres(
    K, rhs, length(rhs)-1;
    tol = tol, maxIter = maxit, M = M, x = x0, out = 1
  )

  (xst, kparams)

end
