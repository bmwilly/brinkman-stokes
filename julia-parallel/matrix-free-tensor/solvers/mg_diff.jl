reload("stokes_flow/agal.jl")
reload("solvers/mg_prolong.jl")
reload("solvers/mg_diff_setup.jl")
reload("solvers/mg_afunbc_diag.jl")
reload("solvers/mg_smooth.jl")

###MG_DIFF GMG preconditioner for diffusion problem
function mg_diff(kparams)

	x = kparams["x"]; y = kparams["y"]
	nc = log2(length(y) - 1)

	# compute new MG data
	h = 2^(1 - nc)

	println("Setting up MG data...")

	# top level
	mgdata = Array{Dict{},int(nc)}
	mgdata[nc] = {
		"name" => "agal",
		"matrix" => (u, params) -> agal(u, params),
		"mat_params" => kparams,
		"prolong" => mg_prolong(2^nc, 2^nc, x, y)
	}

	xn = x; yn = y;
	# loop over remaining levels
	for level = (nc - 1):-1:2
		xn = xn[1:2:end]; yn = yn[1:2:end]
		mparams = mg_diff_setup(xn, yn)
		mgdata[level] = {
			"name" => "mg_afunbc_diag",
			"matrix" => (u, params) -> mg_afunbc_diag(u, params),
			"mat_params" => mparams,
			"prolong" => mg_prolong(2^level, 2^level, xn, yn)
		}
	end

	println("done")

	# MG parameters
	smooth = 1 # point Jacobi
	sweeps = 1; stype = 1

	# pre- and post-smoothing steps
	# npre = 1; npost = 1
	npre = int(input("Number of pre-smoothing steps: "))
	npost = int(input("Number of post-smoothing steps: "))

	# construct smoother
	smooth_data = mg_smooth(mgdata, nc, sweeps, smooth, stype)

	(mgdata, smooth_data, sweeps, stype, npre, npost, nc)

end
