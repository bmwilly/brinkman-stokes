using LinearOperators
using KrylovMethods
using Gadfly
reload("stokes_flow/square_stokes.jl")
reload("stokes_flow/obstacle_stokes.jl")
reload("stokes_flow/brinkman_stokes.jl")
reload("stokes_flow/flowbc.jl")
reload("solvers/mg_diff.jl")
reload("solvers/m_st_mg.jl")

###SOLVE_STOKES solve stokes problem
function solve_stokes(domain)

  @time (if domain == 1
    mats = square_stokes()
  elseif domain == 2
    mats = obstacle_stokes()
  elseif domain == 3
    mats = brinkman_stokes()
  else
    error("invalid domain, please try again")
  end)

  A = mats["A"]; B = mats["B"]; Bx = mats["Bx"]; By = mats["By"];
  f = mats["f"]; g = mats["g"]; xy = mats["xy"]; xyp = mats["xyp"];
  bound = mats["bound"]; x = mats["x"]; y = mats["y"];
  Q = mats["Q"];
  if domain == 3
    P = mats["P"]
  else
    P = zeros(size(A))
  end

  # boundary conditions
  println("imposing (enclosed flow) boundary conditions ...")
  (Ast, Bst, fst, gst) = flowbc(A, B, f, g, xy, bound, domain)
  np = length(gst)
  rhs = vec([fst; gst])

  ## compute solution
  K = [Ast Bst'; Bst spzeros(np, np)]
  M = u -> u

  if domain == 1 || domain == 3
    pc = input("Use GMG preconditioner? (y/n): ")
    if pc == "y"

      nv = size(Ast, 1); np = size(Q, 1); nu = nv/2
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
    end
  end

  tol = 1e-6; maxit = 100;
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
  #     tol = tol, maxIter = maxit, M = M, out = 1
  #   )
  # end
  # )

  if flag == 0
    println("GMRES reached desired tolerance at iteration $(length(resvec))")
  end

  # p = plot(
  #   x = 1:1:length(resvec), y = log(10, resvec), Geom.line,
  #   Guide.xlabel("Iteration"), Guide.ylabel("log_10(residual)")
  # )
  # draw(PNG("graphs/iters.png", 12inch, 6inch), p)

  sol = {
    "K" => K, "Ast" => Ast, "Bst" => Bst, "M" => M, "P" => P,
    "rhs" => rhs, "fst" => fst, "gst" => gst, "xst" => xst
  }

  sol = merge(mats, sol)

end
