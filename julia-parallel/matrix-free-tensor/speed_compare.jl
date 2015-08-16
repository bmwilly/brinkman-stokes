# using Gadfly
# using PyPlot
@everywhere include("hos_homg.jl")
# include("hos_homg.jl")
# @everywhere include("mv_fun.jl")

@everywhere orders = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24];
# orders = [2,3,4,5,6,7,8];
# @everywhere orders = [4,5];
# msizes = [5, 6, 7, 8];

# @everywhere order = 4
@everywhere msize = 5
@everywhere dim = 2

# t = hos_homg(order, msize, dim)

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

# orders
otimes = Float64[];
# # dim = int(input("Dimension: "))
# dim = 2
# if dim == 2
#   msize = 5
# elseif dim == 3
#   msize = 3
# end
for order in orders
  t = hos_homg(order, msize, dim)
  push!(otimes, t)
end
# Gadfly.plot(
#   x = orders, y = otimes, Geom.line,
#   Guide.xlabel("Order"), Guide.ylabel("Time (s)"),
#   Guide.title("100 matvec operations for 32x32 grid for different order elements")
# )

# # mesh sizes
# mtimes = Float64[];
# for msize in msizes
#   t = hos_homg(2, msize-1)
#   push!(mtimes, t)
# end
#
# Gadfly.plot(
#   x = msizes, y = mtimes, Geom.line,
#   Guide.xlabel("Mesh size (log of number of elements)"), Guide.ylabel("Time (s)"),
#   Guide.title("100 matvec operations for Q2 elements on different sized grids")
# )
