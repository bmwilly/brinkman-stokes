reload("stokes_flow/solve_stokes.jl")
reload("stokes_flow/solve_stokes_mg.jl")
reload("graphs/flowplot.jl")
include("helpers/input.jl")

(xst, kparams) = solve_stokes()
println("done")
doplot = input("Create plot? (y/n): ")

if doplot == "y"
  flowplot(xst, kparams)
end
