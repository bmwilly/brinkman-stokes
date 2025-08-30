reload("stokes_flow/agal_diag.jl")

###MG_SMOOTH smoothers for GMG on a square domain
# input
# 	As 				array of dicts containing coefficient matrices
# 	level 			finest level of grid
# 	sweeps			number of directions used for Gauss-Seidel
# 	smooth 			choice of smoother (1:Jacobi, 2:Gauss-Seidel, 3:ILU)
# 	stype 			type of smoother (1:point, 2:line)
# output
# 	Qs 				dict containing smoothing operator in factored form
function mg_smooth(As, level, sweeps, smooth, stype)

	Qs = Array(Dict{}, int(level))

	for i in level:-1:2
		A = As[i]["matrix"]
		params = As[i]["mat_params"]

		# line Gauss-Seidel
		if stype == 2
			#TODO

			# point smoothers
		else

			#ILU
			if smooth == 3
				#TODO

				# point Gauss-Seidel
			elseif smooth == 2
				#TODO

				# point damped Jacobi
			else
				omega = 8 / 9 # relaxation factor
				L1 = (u, params) -> u
				U1 = (u, params) -> (1 / omega) * A(u, params)
				L2 = u -> zeros(length(u), 1)
				U2 = L2
				L3 = L2
				U3 = L2
				L4 = L2
				U4 = L2
			end
		end

		qs = Dict(
			"L1" => L1, "L2" => L2, "L3" => L3, "L4" => L4,
			"U1" => U1, "U2" => U2, "U3" => U3, "U4" => U4,
		)
		Qs[i] = qs

	end # for loop

	return Qs

end
