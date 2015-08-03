using LinearOperators
using KrylovMethods
reload("stokes_flow/square_stokes.jl")
reload("stokes_flow/obstacle_stokes.jl")
reload("stokes_flow/flowbc.jl")
reload("helpers/input.jl")
reload("solvers/mg_diff.jl")
reload("solvers/m_st_mg.jl")

domain = 1

  @time (if domain == 1
    mats = square_stokes()
  elseif domain == 2
    mats = obstacle_stokes()
  else
    error("invalid domain, please try again")
  end)

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
  K = LinearOperator(K)
  M = opEye(length(rhs))

  if domain == 1
    # pc = input("Use GMG preconditioner? (y/n): ")
    pc = "y"
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

      # block MG preconditioner
      M = u -> m_st_mg(u, mparams)
      M = LinearOperator(nv+np, Float64, M)

    end
  end

  # u = linspace(1, size(M, 1), size(M, 1))

  U1 = smooth_data[4]["U1"]
  s1,s2 = size(U1)
  U4u = U1*linspace(1, s1, s1)

  U1 = smooth_data[3]["U1"]
  s1,s2 = size(U1)
  U3u = U1*linspace(1, s1, s1)

  omega = 8/9 # relaxation factor
  Q1 = (1/omega) * spdiagm(diag(A), 0, N, N)

  U1 = smooth_data[2]["U1"]
  s1,s2 = size(U1)
  U2u = U1*linspace(1, s1, s1)


  # Mu = M*u

  # tol = 1e-6;
  # maxit = 100;
  #
  # bnrm2 = norm(rhs);
  # if bnrm2 == 0.0; bnrm2 = 1.0; end
  # err = norm(M*rhs) / bnrm2
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
  #
  # if flag == 0
  #   println("GMRES reached desired tolerance at iteration $(length(resvec))")
  # end
  #
  # sol = merge(mats, {"xst" => xst})
