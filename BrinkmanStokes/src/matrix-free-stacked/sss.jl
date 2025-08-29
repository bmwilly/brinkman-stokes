# using ParallelSparseMatMul  # Package not available in current Julia registry
# Compatibility function for single-threaded execution
using Distributed
share(x) = x
include("stokes_flow/solve_stokes.jl")
include("graphs/flowplot.jl")
include("helpers/helper_functions.jl")

# Compatibility functions for modern Julia
int(x) = Int(x)
input(prompt) = (print(prompt); parse(Int, readline()))

domain = int(input("Choose domain (1/lid-driven cavity, 2/brinkman): "))
msize = int(input("Mesh size: "))
# msize = 1;
# domain = 2;
(sol, kparams) = solve_stokes(domain, msize)
println("done")
doplot = input("Create plot? (y/n): ")

if doplot == "y"
    flowplot(sol, kparams, domain)
end
