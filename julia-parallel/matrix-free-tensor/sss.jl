using ParallelSparseMatMul
include("stokes_flow/solve_stokes.jl")
include("graphs/flowplot.jl")
include("helpers/helper_functions.jl")

domain = int(input("Choose problem (1/lid-driven cavity, 2/brinkman): "))
# domain = 1
(xst, kparams) = solve_stokes(domain)
println("done")
doplot = input("Create plot? (y/n): ")

if doplot == "y"
  flowplot(xst, kparams)
end
