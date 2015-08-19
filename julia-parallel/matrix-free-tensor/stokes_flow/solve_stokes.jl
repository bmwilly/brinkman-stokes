using LinearOperators
using KrylovMethods
reload("stokes_flow/square_stokes.jl")
reload("stokes_flow/brinkman_stokes.jl")
reload("stokes_flow/fbc.jl")
reload("stokes_flow/gbc.jl")
reload("stokes_flow/kfunbc.jl")
reload("helpers/input.jl")
reload("solvers/mg_diff.jl")
reload("solvers/m_st_mg.jl")
reload("stokes_flow/qfun.jl")
reload("stokes_flow/qfun_diag.jl")

###SOLVE_STOKES solve stokes problem
function solve_stokes(domain)

  msize = int(input("Mesh size: "))
  @time (if domain == 1
    kparams = square_stokes(msize)
  elseif domain == 2
    kparams = brinkman_stokes(msize)
  else
    error("invalid domain, please try again")
  end)

  x = kparams["x"]; y = kparams["y"];
  xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"]
  nu = length(xy[:, 1]); nv = 2nu; np = 3length(xyp[:, 1]);

  println("imposing (enclosed flow) boundary conditions ...")
  # get RHS vector
  fst = fbc(kparams)
  gst = gbc(kparams)
  rhs = vec([fst; gst])

  K = u -> kfunbc(u, kparams)
  M = u -> u

  pc = input("Use GMG preconditioner? (y/n): ")
  if pc == "y"

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

    # MG preconditioner
    M = u -> m_st_mg(u, mparams)
  end

  tol = 1e-6; maxit = 100 # gmres parameters
  @time ((xst, flag, err, iter, resvec) = gmres(
    K, rhs, length(rhs);
    tol = tol, maxIter = maxit, M = M, out = 1
  ))
  # bnrm2 = norm(rhs);
  # if bnrm2 == 0.0; bnrm2 = 1.0; end
  # # err = norm(M*rhs) / bnrm2
  # err = norm(M(rhs)) / bnrm2
  # @time(
  # if err < tol
  #   xst = zeros(length(rhs));
  #   flag = NaN; iter = NaN; resvec = NaN
  # else
  #   (xst, flag, err, iter, resvec) = gmres(
  #     K, rhs, length(rhs);
  #     tol = 1e-6, maxIter = maxit, M = M, out = 1
  #   )
  # end
  # )

  if flag == 0
    println("GMRES reached desired tolerance at iteration $(length(resvec))")
  end

  (xst, kparams)

end
