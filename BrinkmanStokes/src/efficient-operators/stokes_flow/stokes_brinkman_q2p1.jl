include("../diffusion/deriv.jl")
include("../diffusion/qderiv.jl")
include("../diffusion/lderiv.jl")

using SparseArrays

###STOKES_Q2P1 Q2-P1 matrix generator
# input
#    xy          Q2 nodal coordinate vector
#    xyp         Q1 nodal coordinate vector
#    mv          Q2 element mapping matrix
# output
#    A           Q2 vector diffusion matrix
#    B           Q2-Q1 divergence matrix
#    Bx          Q2 x-derivative matrix
#    By          Q2 y-derivative matrix
#    f           velocity RHS vector
#    g           pressure RHS vector
function stokes_brinkman_q2p1(grid, K)

	xy = grid["xy"]
	xyp = grid["xyp"]
	mv = grid["mv"]

	nngpt = 9
	x = xy[:, 1]
	y = xy[:, 2]
	xp = xyp[:, 1]
	yp = xyp[:, 2]
	nvtx = length(x)
	nu = 2nvtx
	np = 3length(xp)
	nel = length(mv[:, 1])
	mp = [[1:3:3nel] [2:3:3nel] [3:3:3nel]]

	println("setting up Q2-P1 matrices... ")

	# initialize global matrices
	A = spzeros(nu, nu)
	G = spzeros(nu, nu)
	BBx = spzeros(nvtx, nvtx)
	BBy = spzeros(nvtx, nvtx)
	Bx = spzeros(np, nvtx)
	By = spzeros(np, nvtx)
	B = spzeros(np, nu)
	Q = spzeros(np, np)
	f = zeros(nu, 1)
	g = zeros(np, 1)

	## Gauss point integration rules
	# 3x3 Gauss points
	s = zeros(nngpt, 1)
	t = zeros(nngpt, 1)
	wt = zeros(nngpt, 1)
	gpt = sqrt(0.6)
	s[1] = -gpt
	t[1] = -gpt
	wt[1] = 25 / 81
	s[2] = gpt
	t[2] = -gpt
	wt[2] = 25 / 81
	s[3] = gpt
	t[3] = gpt
	wt[3] = 25 / 81
	s[4] = -gpt
	t[4] = gpt
	wt[4] = 25 / 81
	s[5] = 0.0
	t[5] = -gpt
	wt[5] = 40 / 81
	s[6] = gpt
	t[6] = 0.0
	wt[6] = 40 / 81
	s[7] = 0.0
	t[7] = gpt
	wt[7] = 40 / 81
	s[8] = -gpt
	t[8] = 0.0
	wt[8] = 40 / 81
	s[9] = 0.0
	t[9] = 0.0
	wt[9] = 64 / 81

	# inner loop over elements
	xlv = zeros(nel, 4)
	ylv = zeros(nel, 4)
	for ivtx in 1:4
		xlv[:, ivtx] = x[mv[:, ivtx]]
		ylv[:, ivtx] = y[mv[:, ivtx]]
	end

	# initialize derivative matrices
	ae = zeros(nel, 9, 9)
	bxe = zeros(nel, 3, 9)
	bye = zeros(nel, 3, 9)
	ge = zeros(nel, 9, 9)
	qe = zeros(nel, 3, 3)
	bbxe = zeros(nel, 9, 9)
	bbye = zeros(nel, 9, 9)

	# loop over Gauss points
	for igpt in 1:nngpt
		sigpt = s[igpt]
		tigpt = t[igpt]
		wght = wt[igpt]

		# evaluate derivatives, etc.
		(jac, invjac, phi, dphidx, dphidy) = deriv(sigpt, tigpt, xlv, ylv)
		(psi, dpsidx, dpsidy) = qderiv(sigpt, tigpt, xlv, ylv)
		(chi, dchidx, dchidy) = lderiv(sigpt, tigpt, xlv, ylv)
		# Kigpt = K[:,:,igpt]
		# Kigpt = K;
		# Kigpt = K[1:2:end, 1:2:end];
		for j in 1:9
			for i in 1:9
				ae[:, i, j] += wght * dpsidx[:, i] .* dpsidx[:, j] .* invjac[:]
				ae[:, i, j] += wght * dpsidy[:, i] .* dpsidy[:, j] .* invjac[:]
				# ae[:, i, j] += wght * psi[:, i] .* psi[:, j] .* jac[:] .* Kigpt[:]
				ge[:, i, j] += wght * psi[:, i] .* psi[:, j] .* jac[:]
				bbxe[:, i, j] -= wght * psi[:, i] .* dpsidx[:, j]
				bbye[:, i, j] -= wght * psi[:, i] .* dpsidy[:, j]
			end
			for i in 1:3
				bxe[:, i, j] -= wght * chi[:, i] .* dpsidx[:, j]
				bye[:, i, j] -= wght * chi[:, i] .* dpsidy[:, j]
			end
		end

		for j in 1:3
			for i in 1:3
				qe[:, i, j] += wght * chi[:, i] .* chi[:, j] .* jac[:]
			end
		end
	end # end of Gauss point loop

	## element assembly into global matrices
	# component velocity matrices
	for krow in 1:9
		nrow = mv[:, krow]
		for kcol in 1:9
			ncol = mv[:, kcol]
			A += sparse(nrow, ncol, ae[:, krow, kcol], nu, nu)
			A += sparse(nrow .+ nvtx, ncol .+ nvtx, ae[:, krow, kcol], nu, nu)
			G += sparse(nrow, ncol, ge[:, krow, kcol], nu, nu)
			G += sparse(nrow .+ nvtx, ncol .+ nvtx, ge[:, krow, kcol], nu, nu)
			BBx += sparse(nrow, ncol, bbxe[:, krow, kcol], nvtx, nvtx)
			BBy += sparse(nrow, ncol, bbye[:, krow, kcol], nvtx, nvtx)
		end
		for kcol in 1:3
			ncol = mp[:, kcol]
			Bx += sparse(ncol, nrow, bxe[:, kcol, krow], np, nvtx)
			By += sparse(ncol, nrow, bye[:, kcol, krow], np, nvtx)
		end
	end

	# vector velocity matrices
	B = [Bx By]

	# pressure matrices
	for krow in 1:3
		nrow = mp[:, krow]
		for kcol in 1:3
			ncol = mp[:, kcol]
			Q += sparse(nrow, ncol, qe[:, krow, kcol], np, np)
		end
	end

	println("done")

	A += K .* G # brinkman
	# G = K .* G

	mats = Dict(
		"A" => A, # velocity matrix
		"B" => B, # vector velocity matrix
		"G" => G, # mass matrix
		"Q" => Q, # pressure matrix
		"Bx" => BBx,
		"By" => BBy,
		"P" => K, # permeability tensor
		"f" => f,
		"g" => g,
	)
	return mats

end
