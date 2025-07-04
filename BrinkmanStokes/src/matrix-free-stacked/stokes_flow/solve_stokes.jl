using KrylovMethods
using DelimitedFiles
using PyPlot
include("../../output_utils.jl")
include("square_stokes.jl")
include("brinkman_stokes.jl")
include("fbc.jl")
include("gbc.jl")
include("kfunbc.jl")
include("qfun.jl")
include("qfun_diag.jl")
include("../solvers/mg_diff.jl")
include("../solvers/m_st_mg.jl")

###SOLVE_STOKES solve stokes problem
function solve_stokes(domain, msize)

	@time (
		if domain == 1
			kparams = square_stokes(msize)
		elseif domain == 2
			kparams = brinkman_stokes(msize)
		else
			error("invalid domain, please try again")
		end
	)

	x = kparams["x"]
	y = kparams["y"]
	xy = kparams["xy"]
	xyp = kparams["xyp"]
	mv = kparams["mv"]
	bound = kparams["bound"]
	nu = length(xy[:, 1])
	nv = 2nu
	np = 3length(xyp[:, 1])

	println("imposing (enclosed flow) boundary conditions ...")
	# get RHS vector
	fst = fbc(domain, kparams)
	gst = gbc(domain, kparams)
	rhs = vec([fst; gst])

	K = u -> kfunbc(u, kparams, domain)
	M = u -> u

	pc = input("Use GMG preconditioner? (y/n): ")
	if pc == "y"

		nu = length(xy[:, 1])
		nv = 2nu
		np = 3length(xyp[:, 1])

		# set up GMG preconditioner
		(mgdata, smooth_data, sweeps, stype, npre, npost, nc) = mg_diff(kparams)

		mparams = {
			"nv" => nv,
			"np" => np,
			"Q" => u -> qfun(u, kparams),
			"Qdiag" => u -> qfun_diag(u, kparams),
			"mgdata" => mgdata,
			"smooth_data" => smooth_data,
			"nc" => nc,
			"npre" => npre,
			"npost" => npost,
			"sweeps" => sweeps,
		}

		# MG preconditioner
		M = u -> m_st_mg(u, mparams)
	end

	restrt = min(5000, length(rhs))
	tol = 1e-6
	maxIter = 100
	@time ((xst, flag, err, iter, resvec) = gmres(
		K, rhs, restrt;
		tol = tol, maxIter = maxIter, M = M, out = 1,
	))

	if flag == 0
		println("GMRES reached desired tolerance at iteration $(length(resvec))")
	end

	# Create convergence plot using unified output system
	convergence_file = get_output_file("matrix-free-stacked", domain, msize, "convergence.png"; subdir = "plots")
	figure(figsize = (10, 6))
	plot(1:length(resvec), log10.(resvec), "b-", linewidth = 2)
	xlabel("Iteration")
	ylabel("log₁₀(residual)")
	title("GMRES Convergence (Domain=$domain, Size=$msize)")
	grid(true)
	savefig(convergence_file, dpi = 150, bbox_inches = "tight")
	println("Convergence plot saved to: $(convergence_file)")
	close("all")

	# Write CSV file using unified output system
	outfile = get_output_file("matrix-free-stacked", domain, msize, "convergence_data.csv")
	DelimitedFiles.writedlm(outfile, log.(10, resvec))

	# Add msize to kparams for plotting
	kparams["msize"] = msize

	(xst, kparams)

end
