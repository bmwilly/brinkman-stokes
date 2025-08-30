include("afunbc.jl")
include("ho_afun.jl")
include("bfunbc.jl")
include("btfunbc.jl")

###KFUNBC matrix-free stiffness operator
function kfunbc(u, kparams)
	if domain == 1
		af(u, kparams) = afunbc(u, kparams)
	elseif domain == 2
		af(u, kparams) = ho_afun(u, kparams)
	end

	u1 = copy(u)
	u2 = copy(u)
	u3 = copy(u)
	w = af(u1, kparams) + bfunbc(u2, kparams) + btfunbc(u3, kparams)
	return w
end
