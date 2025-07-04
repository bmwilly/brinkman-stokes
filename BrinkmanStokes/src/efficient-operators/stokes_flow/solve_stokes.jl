using DelimitedFiles
using KrylovMethods
using PyPlot
# using IterativeSolvers
include("square_stokes.jl")
include("brinkman_stokes.jl")
include("flowbc.jl")
include("../solvers/mg_diff.jl")
include("../solvers/m_st_mg.jl")

###SOLVE_STOKES solve stokes problem
function solve_stokes(domain::Int, msize::Int)

	@time (
		if domain == 1
			mats = square_stokes(msize)
		elseif domain == 2
			mats = brinkman_stokes(msize)
		else
			error("invalid domain, please try again")
		end
	)

	A = mats["A"]
	B = mats["B"]
	Bx = mats["Bx"]
	By = mats["By"]
	f = mats["f"]
	g = mats["g"]
	xy = mats["xy"]
	xyp = mats["xyp"]
	bound = mats["bound"]
	x = mats["x"]
	y = mats["y"]
	Q = mats["Q"]
	msize = mats["msize"]
	if domain == 2
		kappa = mats["kappa"]
	else
		kappa = zeros(length(x) * length(y))
	end

	# boundary conditions
	println("imposing (enclosed flow) boundary conditions ...")
	(Ast, Bst, fst, gst) = flowbc(A, B, f, g, xy, bound, domain)
	np = length(gst)
	rhs = vec([fst; gst])

	# compute solution
	K = [Ast Bst'; Bst spzeros(np, np)]
	M = u -> u

	# TODO:
	# pc = user_input("Use GMG preconditioner? (y/n): ")
	# pc = "y"
	pc = "n"
	if pc == "y"

		nv = size(Ast, 1)
		np = size(Q, 1)
		nu = Int(nv / 2)
		Agal = Ast[1:nu, 1:nu]
		(mgdata, smooth_data, sweeps, stype, npre, npost, nc) = mg_diff(x, y, Agal)

		mparams = Dict(
			"nv" => nv,
			"Q" => Q,
			"mgdata" => mgdata,
			"smooth_data" => smooth_data,
			"nc" => nc,
			"npre" => npre,
			"npost" => npost,
			"sweeps" => sweeps,
		)

		# block GMG preconditioner
		M = u -> m_st_mg(u, mparams)
	end

	restrt = min(5000, length(rhs))
	tol = 1e-6
	maxIter = 100
	@time ((xst, flag, err, iter, resvec) = gmres(
		K, rhs, restrt;
		tol = tol, maxIter = maxIter, M = M, out = 1,
	))

	# m_st_mg!(u, unneeded, mparams) = m_st_mg(u, mparams)
	# M = MatrixFcn{Float64}(size(K,1), size(K,2), (u, unneeded) -> m_st_mg!(u, unneeded, mparams))
	# @time ((xst, convHist) = IterativeSolvers.gmres(K, rhs, M; tol = tol, restart = restrt))
	# xst,convHist = gmres(K, rhs)

	if flag == 0
		println("GMRES reached desired tolerance at iteration $(length(resvec))")
	end

	# Define repo root for output files and plots
	repo_root = joinpath(@__DIR__, "..")

	# Create convergence plot
	plots_dir = joinpath(repo_root, "output", "plots")
	mkpath(plots_dir)

	figure(figsize = (10, 6))
	plot(1:length(resvec), log10.(resvec), "b-", linewidth = 2)
	xlabel("Iteration")
	ylabel("log₁₀(residual)")
	title("GMRES Convergence")
	grid(true)
	convergence_file = joinpath(plots_dir, "convergence_msize$(msize).png")
	savefig(convergence_file, dpi = 150, bbox_inches = "tight")
	println("Convergence plot saved to: $(convergence_file)")
	close("all")

	# Write to output folder in repo root, regardless of current working directory
	output_dir = joinpath(repo_root, "output")
	mkpath(output_dir)  # Create output directory if it doesn't exist
	outfile = joinpath(output_dir, "brinkman_iters$(msize).csv")
	DelimitedFiles.writedlm(outfile, log.(10, resvec))

	sol = Dict(
		"K" => K, "Ast" => Ast, "Bst" => Bst, "M" => M, "kappa" => kappa,
		"rhs" => rhs, "fst" => fst, "gst" => gst, "xst" => xst,
		"resvec" => resvec, "err" => err, "iter" => iter,
	)

	sol = merge(mats, sol)
end
