reload("stokes_flow/solve_stokes.jl")
reload("stokes_flow/solve_stokes_mg.jl")
reload("graphs/flowplot.jl")
include("helpers/input.jl")

domain = int(input("Choose problem (1/lid-driven cavity, 2/brinkman): "))
(xst, kparams) = solve_stokes(domain)
println("done")
doplot = input("Create plot? (y/n): ")

if doplot == "y"
  flowplot(xst, kparams)
end
