using ParallelSparseMatMul
include("stokes_flow/solve_stokes.jl")
include("graphs/flowplot.jl")
include("helpers/helper_functions.jl")

domain = int(input("Choose domain (1/lid-driven cavity, 2/brinkman): "))
msize = int(input("Mesh size: "))
# msize = 1;
# domain = 2;
(sol, kparams) = solve_stokes(domain, msize)
println("done")
doplot = input("Create plot? (y/n): ")

if doplot == "y"
  flowplot(sol, kparams)
end
