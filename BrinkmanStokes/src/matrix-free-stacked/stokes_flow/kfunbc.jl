using LinearOperators
include("afunbc.jl")
include("ho_afun.jl")
include("bfunbc.jl")
include("btfunbc.jl")

###KFUNBC matrix-free stiffness operator
function kfunbc(u, kparams, domain)

	xy = kparams["xy"]
	xyp = kparams["xyp"]
	nvtx = length(xy[:, 1])
	nu = 2nvtx
	np = 3length(xyp[:, 1])

	if domain == 1
		A = LinearOperator(nu + np, Float64, u -> afunbc(u, kparams))
	elseif domain == 2
		A = LinearOperator(nu + np, Float64, u -> ho_afun(u, kparams))
	end
	B = LinearOperator(nu + np, Float64, u -> bfunbc(u, kparams))
	Bt = LinearOperator(nu + np, Float64, u -> btfunbc(u, kparams))

	return A * u + B * u + Bt * u

end
