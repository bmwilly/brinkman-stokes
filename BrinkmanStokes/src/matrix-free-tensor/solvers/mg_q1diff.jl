using SparseArrays

###MG_Q1DIFF bilinear diffusion matrix generator for GMG
# input
# 	xy 				vertex coordinate vector
# 	ev 				element mapping matrix
# output
# 	a 				stiffness matrix
# 	r 				dummy variable
# 	f 				rhs vector
function mg_q1diff(xy, ev)

    x = xy[:, 1]
    y = xy[:, 2]
    nvtx = length(x)
    nel = length(ev[:, 1])
    lx = maximum(x) - minimum(x)
    ly = maximum(y) - minimum(y)
    hx = maximum(diff(x))
    hy = maximum(diff(y))

    # initialize global matrices
    # a = spzeros(nvtx, nvtx)
    # r = spzeros(nvtx, nvtx)
    f = zeros(nvtx, 1)
    ae = zeros(nel, 4, 4)

    # set up 2x2 Gauss points
    s = zeros(4, 1)
    t = zeros(4, 1)
    gpt = 1.0 / sqrt(3.0)
    s[1] = -gpt
    t[1] = -gpt
    s[2] = gpt
    t[2] = -gpt
    s[3] = gpt
    t[3] = gpt
    s[4] = -gpt
    t[4] = gpt

    # inner loop over elements
    xlv = zeros(nel, 4)
    ylv = zeros(nel, 4)
    for ivtx in 1:4
        xlv[:, ivtx] = x[ev[:, ivtx]]
        ylv[:, ivtx] = y[ev[:, ivtx]]
    end

    # loop over 2x2 Gauss points
    for igpt in 1:4
        sigpt = s[igpt]
        tigpt = t[igpt]

        # evaluate derivatives, etc
        (jac, invjac, phi, dphidx, dphidy) = deriv(sigpt, tigpt, xlv, ylv)
        for j in 1:4
            for i in 1:4
                ae[:, i, j] += dphidx[:, i] .* dphidx[:, j] .* invjac[:]
                ae[:, i, j] += dphidy[:, i] .* dphidy[:, j] .* invjac[:]
            end
        end
    end # Gauss point loop

    return ae, f

    # # assemble global matrix and source vector
    # for krow = 1:4
    # 	nrow = int(ev[:, krow])
    # 	for kcol = 1:4
    # 		ncol = int(ev[:, kcol])
    # 		a += sparse(nrow, ncol, ae[:, krow, kcol], nvtx, nvtx)
    # 	end
    # end
    #
    # (a, r, f)

end

# ###MG_AFUN constructs matrix-free bilinear diffusion operator for GMG
# function mg_afun(u, mparams)
#
# 	ev = mparams["ev"]
# 	ae = mparams["ae"]
# 	aes = squeeze(ae[1, :, :], 1)
#
# 	for e = 1:nel
# 		ind = ev[e, :]'
# 		we = aes * u[ind]
# 		w[ind] += we
# 	end
#
# 	vec(w)
#
# end
