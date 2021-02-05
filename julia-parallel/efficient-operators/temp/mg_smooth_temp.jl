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

nv = size(Ast, 1); np = size(Q, 1); nu = nv / 2
Agal = Ast[1:nu, 1:nu]

nc = log2(length(y) - 1)

  # compute new MG data
h = 2^(1 - nc)

println("Setting up MG data...")

  # top level
mgdata = Array{Dict{},int(nc)}
mgdata[nc] = {
    "matrix" => Agal,
    "prolong" => mg_prolong(2^nc, 2^nc, x, y)
  }

xn = x; yn = y
  # loop over remaining levels
for level = (nc - 1):-2:2
    xn = xn[1:2:end]; yn = yn[1:2:end]
    println("size(xn) $(size(xn))")
    mgdata[level] = {
      "matrix" => mg_diff_setup(xn, yn),
      "prolong" => mg_prolong(2^level, 2^level, xn, yn)
    }
end

println("done")

  # MG parameters
smooth = 1 # point Jacobi
sweeps = 1; stype = 1

  # pre- and post-smoothing steps
npre = 1; npost = 1






  ## mg_smooth
Qs = Array(Dict{}, int(level))

