reload("solvers/mg_pre.jl")
reload("solvers/mg_post.jl")
reload("helpers/gmres2.jl")

###MG_ITER performs one GMG iteration
# input
# 	As				coefficient matrix
# 	x0 				initial iterate
# 	f 				rhs
# 	smooth_data		dict containing smoothing operators
# 	level			grid level
# 	npre 			number of presmoothing steps
# 	npost 			number of postsmoothing steps
# 	sweeps			type of smoothing strategy used for Gauss-Seidel smoothing
# output
# 	x 				result of multigrid steps
function mg_iter(As, x0, f, smooth_data, level, npre, npost, sweeps)

	A = As[level]["matrix"];
	params = As[level]["mat_params"]
	name = As[level]["name"]
	P = As[level]["prolong"]

	if level == 2
		x,ignore = gmres(u -> A(u, params), f, length(f))
		# x,flag,err,iter,resvec = gmres(u -> A(u, params), f, length(f))
	else
		# presmooth
		x = mg_pre(A, params, x0, f, npre, smooth_data, level, sweeps)

		# restrict residual
		r = f - A(x, params)
		rc = P'*r

		# coarse grid correction
		cc = mg_iter(As, zeros(size(rc)), rc, smooth_data, level-1, npre, npost, sweeps)
		x += P*cc

		# postsmooth
		x = mg_post(A, params, x, f, npost, smooth_data, level, sweeps)
	end

	x
end
