using BrinkmanStokes

include("efficient-operators/stokes_flow/solve_stokes.jl")
include("efficient-operators/graphs/flowplot.jl")
include("efficient-operators/helpers/helper_functions.jl")

domain = user_input("Choose domain (1/lid-driven cavity, 2/brinkman): ")
msize = user_input("Mesh size: ")

# Add option for multigrid preconditioning for large problems
use_mg = false
if msize >= 5
	mg_choice = user_input("Use multigrid preconditioning for better performance? (y/n): ")
	use_mg = (mg_choice == "y" || mg_choice == "Y")
end

sol = solve_stokes(domain, msize; use_mg = use_mg)

println("done")
flowplot(sol, domain)
