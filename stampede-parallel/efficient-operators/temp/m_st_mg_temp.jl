reload("solvers/gmres.jl")
reload("solvers/mg_diff.jl")
reload("solvers/m_st_mg.jl")

###SOLVE_STOKES solve stokes problem
# function solve_stokes_mg()

  (A, B, G, Q, Bx, By, f, g, x, y, xy, xyp, bound) = square_stokes()

  # boundary conditions
  println("imposing (enclosed flow) boundary conditions ...")
  (Ast, Bst, fst, gst) = flowbc(A, B, f, g, xy, bound)
  np = length(gst)

  ## STOKES problem

  # set parameters
  tol = 1e-6; 
  maxit = 100;

  # specify coefficient matrix
  K = [Ast Bst'; Bst spzeros(np, np)]

  # construct preconditioner
  println("multigrid preconditioning...")

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


  nv = mparams["nv"]; np = size(mparams["Q"], 1); nu = nv/2
  x_it = rand(nv + np, 1)
  rvx = x_it[1:nu]; rvy = x_it[nu+1:nv]; rp = x_it[nv+1:nv+np]

  # zvx = mg_iter(
  #   mparams["mgdata"], zeros(size(rvx)), rvx, 
  #   mparams["smooth_data"], mparams["nc"], 
  #   mparams["npre"], mparams["npost"], mparams["sweeps"]
  # )

As = mparams["mgdata"]; x0 = zeros(size(rvx)); f = rvx; 
smooth_data = mparams["smooth_data"]; level = mparams["nc"]; 
npre = mparams["npre"]; npost = mparams["npost"]; sweeps = mparams["sweeps"];

A = As[level]["matrix"]; P = As[level]["prolong"]; 

# x = mg_pre(A, x0, f, npre, smooth_data, level, sweeps)

xs0 = x0; ns = npre; Qs = smooth_data; 

L1 = Qs[level]["L1"]

  # M = u -> m_st_mg(u, mparams)

  # println("GMRES iteration...")

  # # zero initial guess
  # x0 = vec(zeros(size([fst; gst])))

  # tic()
  # (xst, flag, err, iter, resvec, cnt) = gmres(
  #   K, vec([fst; gst]), length([fst; gst]); 
  #   tol = tol, maxIter = maxit, M = M, x = x0, out = 1
  # )
  # etoc = toc()
  # println("gmres converged at iteration $(cnt) to a solution with relative residual $(err)")
  # println("Stokes system solved in $(etoc) seconds")

  # (xst, By, Bx, A, xy, xyp, bound)

# end