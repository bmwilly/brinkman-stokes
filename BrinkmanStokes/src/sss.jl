using BrinkmanStokes

include("efficient-operators/stokes_flow/solve_stokes.jl")
include("efficient-operators/graphs/flowplot.jl")
include("efficient-operators/helpers/helper_functions.jl")

domain = user_input("Choose domain (1/lid-driven cavity, 2/brinkman): ")
msize = user_input("Mesh size: ")
sol = solve_stokes(domain, msize)

println("done")
flowplot(sol, domain)
