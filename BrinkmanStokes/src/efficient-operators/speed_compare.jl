push!(LOAD_PATH, "$(homedir())/Documents/brinkman-stokes/julia-parallel/efficient-operators")
include("hos_homg.jl")
include("mv_fun.jl")

# orders = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24];
# orders = [15,16,17,18,19,20];
# orders = [21,22,23,24];
# orders = [2,3,4,5,6,7,8,9,10];
# orders = [11,12,13,14,15,16,17,18,19,20,21,22,23,24];
# msizes = [5,6,7];
# msizes = [9]

# msize = 9
msize = user_input("Mesh size: ")
mats = square_stokes(msize)
A = mats["A"]; B = mats["B"]; Bx = mats["Bx"]; By = mats["By"];
f = mats["f"]; g = mats["g"]; xy = mats["xy"]; xyp = mats["xyp"];
bound = mats["bound"]; x = mats["x"]; y = mats["y"];
Q = mats["Q"]

# boundary conditions
(Ast, Bst, fst, gst) = flowbc(A, B, f, g, xy, bound, 1)
np = length(gst)
rhs = vec([fst; gst])

## compute solution
K = [Ast Bst'; Bst spzeros(np, np)]
n, m = size(K)

tic(); for cnt in 1:100
    u = rand(size(Ast, 1)); w = Ast * u
end
etoc = toc()

# # mesh sizes
# mtimes = Float64[]
# for msize in msizes
#   t = mv_fun(msize)
#   push!(mtimes, t)
# end
# Gadfly.plot(
#   x = msizes, y = mtimes, Geom.line,
#   Guide.xlabel("Mesh size (log of number of elements)"), Guide.ylabel("Time (s)"),
#   Guide.title("100 matvec operations for Q2 elements on different sized grids")
# )

# otimes = Float64[];
# dim = user_input("Dimension: ")
# if dim == 2
#   msize = 5
# elseif dim == 3
#   msize = 3
# end
# for order in orders
#   t = hos_homg(order, msize, dim)
#   push!(otimes, t)
# end
# Gadfly.plot(
#   x = orders, y = otimes, Geom.line,
#   Guide.xlabel("Order"), Guide.ylabel("Time (s)"),
#   Guide.title("100 matvec operations for 32x32 grid for different order elements")
# )

# order = 2;
# dim = 2;
# mtimes = Float64[];
# for msize in msizes
#   t = hos_homg(order, msize, dim)
#   push!(mtimes, t)
# end
# Gadfly.plot(
#   x = msizes, y = mtimes, Geom.line,
#   Guide.xlabel("Mesh size (log of number of elements)"), Guide.ylabel("Time (s)"),
#   Guide.title("100 matvec operations for Q2 elements on different sized grids")
# )
