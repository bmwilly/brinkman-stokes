push!(LOAD_PATH, "$(homedir())/Documents/brinkman-stokes/julia-parallel/matrix-free-tensor")
# using Gadfly
# using PyPlot
include("hos_homg.jl")
include("mv_fun.jl")

# p = int(input("Total number of processors: "))
# addprocs(p-1)

# orders = [21,22,23,24];
# msizes = [5, 6, 7, 8];
msizes = [2,3,4,5,6,7,8];
# msizes = [8]

# msize = 9
# msize = int(input("Mesh size: "))
# kparams = square_stokes(msize)
# xy = kparams["xy"]; xyp = kparams["xyp"]
# nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])

# tic(); for cnt = 1:100; u = Base.shmem_rand(nu+np); w = afunbc(u, kparams); end; etoc = toc()
# tic(); for cnt = 1:100; u = Base.shmem_rand(nu+np); w = afunbc(u, kparams); end; etoc = toc()


# mesh sizes
mtimes = Float64[]
for msize in msizes
  t = mv_fun(msize)
  push!(mtimes, t)
end
@show mtimes
# Gadfly.plot(
#   x = msizes, y = mtimes, Geom.line,
#   Guide.xlabel("Mesh size (log of number of elements)"), Guide.ylabel("Time (s)"),
#   Guide.title("100 matvec operations for Q2 elements on different sized grids")
# )

# orders
# otimes = Float64[];
# dim = int(input("Dimension: "))
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
