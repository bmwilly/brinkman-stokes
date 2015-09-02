push!(LOAD_PATH, "$(homedir())/Documents/brinkman-stokes/julia-parallel/matrix-free-stacked")
# using Gadfly
# using PyPlot
include("hos_homg.jl")
include("mv_fun.jl")

# orders = [2, 3, 4];
msize = [1]
# msizes = [2,3,4,5,6]
# msizes = [2,3,4]

# mesh sizes
# mtimes = Float64[]
# for msize in msizes
#   t = mv_fun(msize)
#   push!(mtimes, t)
# end
# @show mtimes
# Gadfly.plot(
#   x = msizes, y = mtimes, Geom.line,
#   Guide.xlabel("Mesh size (log of number of elements)"), Guide.ylabel("Time (s)"),
#   Guide.title("100 matvec operations for Q2 elements on different sized grids")
# )

# orders
# otimes = Float64[];
# for order in orders
#   t = hos_homg(order, )
#   push!(otimes, t)
# end
#
# Gadfly.plot(
#   x = orders, y = otimes, Geom.line,
#   Guide.xlabel("Order"), Guide.ylabel("Time (s)"),
#   Guide.title("100 matvec operations for 32x32 grid for different order elements")
# )

# mesh sizes
order = 2; dim = 3;
mtimes = Float64[];
for msize in msizes
  t = hos_homg(order, msize, dim)
  push!(mtimes, t)
end
#
# Gadfly.plot(
#   x = msizes, y = mtimes, Geom.line,
#   Guide.xlabel("Mesh size (log of number of elements)"), Guide.ylabel("Time (s)"),
#   Guide.title("100 matvec operations for Q2 elements on different sized grids")
# )
