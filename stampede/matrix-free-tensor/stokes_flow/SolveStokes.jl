
module SolveStokes
# using KrylovMethods
# include("square_stokes.jl")
# include("fbc.jl")
# include("gbc.jl")
# include("../helpers/gmres.jl")
# require("stokes_flow/square_stokes.jl")
# require("stokes_flow/fbc.jl")
# require("stokes_flow/gbc.jl")
# require("helpers/gmres.jl")
reload("stokes_flow/square_stokes.jl")
reload("stokes_flow/fbc.jl")
reload("stokes_flow/gbc.jl")
reload("helpers/gmres.jl")

export solve_stokes 

###SOLVE_STOKES solve stokes problem
function solve_stokes()

  tic()
  (xy, xyp, mv, bound, ae, bxe, bye, ge, qe, bbxe, bbye) = square_stokes()
  toc1 = toc()
  println("square stokes: $(toc1)")

  x = xy[:, 1]
  y = xy[:, 2]
  xp = xyp[:, 1]
  yp = xyp[:, 2]
  nvtx = length(x)
  nu = 2nvtx
  np = 3length(xp)
  nel = length(mv[:, 1])
  mp = [[1:3:3nel]' [2:3:3nel]' [3:3:3nel]']

  println("imposing (enclosed flow) boundary conditions ...")

  ## compute solution
  # get RHS vector
  tic()
  fst = fbc(xy, xyp, mv, bound, ae, bxe, bye)
  gst = gbc(xy, xyp, mv, bound, bxe, bye)
  toc2 = toc()
  println("fbc gbc: $(toc2)")

  tic()
  (xst, flag, err, iter, resvec, cnt) = gmres(
    u -> kfunbc(u, xy, xyp, mv, bound, ae, bxe, bye), [fst; gst], 1000; 
    tol = 1e-6, maxIter = 300, M = u -> mfun(u, xy, xyp, mv, bound, ge, qe), out = 1
  )
  # pgmres =
  # (xst, flag, err, iter, resvec, cnt) = pmap(gmres, )
  etoc = toc()
  println("gmres converged at iteration $(cnt) to a solution with relative residual $(err)")
  println("Stokes system solved in $(etoc) seconds")
  (xst, xy, xyp, mv, bound, ae, bxe, bye, bbxe, bbye)

end

end 