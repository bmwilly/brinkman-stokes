reload("stokes_flow/solve_stokes.jl")
reload("graphs/flowplot.jl")
reload("helpers/helper_functions.jl")

domain = int(input("Choose domain (1/lid-driven cavity, 2/channel with obstacles, 3/brinkman): "))
sol = solve_stokes(domain)
println("done")
flowplot(sol, domain)
