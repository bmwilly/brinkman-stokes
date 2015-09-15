using KrylovMethods
# using IterativeSolvers
include("square_stokes.jl")
include("brinkman_stokes.jl")
include("flowbc.jl")
include("../solvers/mg_diff.jl")
include("../solvers/m_st_mg.jl")

###SOLVE_STOKES solve stokes problem
function solve_stokes(domain, msize)

  @time (if domain == 1
    mats = square_stokes(msize)
  elseif domain == 2
    mats = brinkman_stokes(msize)
  else
    error("invalid domain, please try again")
  end)

  A = mats["A"]; B = mats["B"]; Bx = mats["Bx"]; By = mats["By"];
  f = mats["f"]; g = mats["g"]; xy = mats["xy"]; xyp = mats["xyp"];
  bound = mats["bound"]; x = mats["x"]; y = mats["y"];
  Q = mats["Q"]; msize = mats["msize"];
  if domain == 2; kappa = mats["kappa"]; else; kappa = zeros(length(x)*length(y)); end;

  # boundary conditions
  println("imposing (enclosed flow) boundary conditions ...")
  (Ast, Bst, fst, gst) = flowbc(A, B, f, g, xy, bound, domain)
  np = length(gst)
  rhs = vec([fst; gst])

  ## compute solution
  K = [Ast Bst'; Bst spzeros(np, np)]
  M = u -> u

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

  restrt = min(5000, length(rhs)); tol = 1e-6; maxIter = 100
  @time ((xst, flag, err, iter, resvec) = gmres(
    K, rhs, restrt;
    tol = tol, maxIter = maxIter, M = M, out = 1
  ))

  # m_st_mg!(u, unneeded, mparams) = m_st_mg(u, mparams)
  # M = MatrixFcn{Float64}(size(K,1), size(K,2), (u, unneeded) -> m_st_mg!(u, unneeded, mparams))
  # @time ((xst, convHist) = IterativeSolvers.gmres(K, rhs, M; tol = tol, restart = restrt))
  # xst,convHist = gmres(K, rhs)

  if flag == 0
    println("GMRES reached desired tolerance at iteration $(length(resvec))")
  end

  # p = Gadfly.plot(
  #   x = 1:1:length(resvec), y = log(10, resvec), Geom.line,
  #   Guide.xlabel("Iteration"), Guide.ylabel("log_10(residual)")
  # )
  # outfile = string("graphs/brinkman_iters", msize, ".png")
  # Gadfly.draw(PNG(outfile, 12inch, 6inch), p)
  #
  outfile = string("output/brinkman_iters", msize, ".csv")
  writecsv(outfile, log(10, resvec))

  sol = {
    "K" => K, "Ast" => Ast, "Bst" => Bst, "M" => M, "kappa" => kappa,
    "rhs" => rhs, "fst" => fst, "gst" => gst, "xst" => xst,
    "resvec" => resvec, "err" => err, "iter" => iter
  }

  sol = merge(mats, sol)
end
