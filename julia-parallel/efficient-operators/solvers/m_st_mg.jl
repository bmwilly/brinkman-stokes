reload("solvers/mg_iter.jl")

###M_ST_MG block MG preconditioner for Stokes equations
# input
# 	x_it 			operand for preconditioning operator
# 	mparams 		structure defining MG preconditioning operator
# output
# 	x_it 			result of MG preconditioning operation
function m_st_mg(x_it, mparams)

	nv = mparams["nv"]; np = size(mparams["Q"], 1); nu = nv/2
	rvx = x_it[1:nu];
	rvy = x_it[nu+1:nv];
	rp = x_it[nv+1:nv+np]

	zvx = mg_iter(
		mparams["mgdata"], zeros(size(rvx)), rvx,
		mparams["smooth_data"], mparams["nc"],
		mparams["npre"], mparams["npost"], mparams["sweeps"]
	)

	zvy = mg_iter(
		mparams["mgdata"], zeros(size(rvy)), rvy,
		mparams["smooth_data"], mparams["nc"],
		mparams["npre"], mparams["npost"], mparams["sweeps"]
	)

	zp = diag(mparams["Q"]) .\ rp

	x_it = [zvx; zvy; zp]

end
