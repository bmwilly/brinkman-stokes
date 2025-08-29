include("mg_pre.jl")
include("mg_post.jl")

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

    A = As[level]["matrix"]
    P = As[level]["prolong"]

    if level == 2
        x = A \ f
    else
        # presmooth
        x = mg_pre(A, x0, f, npre, smooth_data, level, sweeps)

        # restrict residual
        r = f - A * x
        rc = P' * r

        # coarse grid correction
        cc = mg_iter(As, zeros(size(rc)), rc, smooth_data, level - 1, npre, npost, sweeps)
        x = x + P * cc

        # postsmooth
        x = mg_post(A, x, f, npost, smooth_data, level, sweeps)
    end

    return x

end
