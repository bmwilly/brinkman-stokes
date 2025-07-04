###MG_AFUN constructs matrix-free bilinear diffusion operator for GMG
function mg_afun(u, mparams)

	ev = mparams["ev"]; ae = mparams["ae"]
	nel = length(ev[:, 1])
	aes = squeeze(ae[1, :, :], 1)
	w = zeros(length(u), 1)

	for e = 1:nel
		ind = ev[e, :]'
		# println(ind)
		we = aes * u[ind]
		w[ind] += we
	end

	vec(w)

end
