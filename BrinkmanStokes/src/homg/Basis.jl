module Basis
export polynomial, gradient, gauss, gll

function polynomial(x::Array, alpha::Real, beta::Real, N::Integer)
	# function [P] = basis.polynomial(x,alpha,beta,N)
	# Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
	# (alpha+beta <> -1) at points x for order N and returns P[1:length(xp))]
	# Note : They are normalized to be orthonormal.
	# Turn points into row if needed.
	xp = reshape(x, size(x)[1], 1)
	dims = size(xp)
	if dims[2] == 1
		xp = xp'
	end
	PL = zeros(N + 1, length(xp))
	# Initial values P_0(x) and P_1(x)
	gamma0 = 2^(alpha + beta + 1) / (alpha + beta + 1) * gamma(alpha + 1) * gamma(beta + 1) / gamma(alpha + beta + 1)
	PL[1,:] = 1.0 / sqrt(gamma0)
	if N == 0
		P = PL'
		return P
	end
	gamma1 = (alpha + 1) * (beta + 1) / (alpha + beta + 3) * gamma0
	PL[2,:] = ((alpha + beta + 2) * xp / 2 + (alpha - beta) / 2) / sqrt(gamma1)
	if N == 1
		P = PL[N + 1,:]'
		return P
	end
	# Repeat value in recurrence.
	aold = 2 / (2 + alpha + beta) * sqrt((alpha + 1) * (beta + 1) / (alpha + beta + 3))
	# Forward recurrence using the symmetry of the recurrence.
	for i = 1:N - 1
		h1 = 2 * i + alpha + beta
		anew = 2 / (h1 + 2) * sqrt((i + 1) * (i + 1 + alpha + beta) * (i + 1 + alpha) * (i + 1 + beta) / (h1 + 1) / (h1 + 3))
		bnew = - (alpha^2 - beta^2) / h1 / (h1 + 2)
		PL[i + 2,:] = 1 / anew * ( -aold * PL[i,:] + (xp - bnew) .* PL[i + 1,:])
		aold = anew
	end
	P = PL[N + 1,:]'
	return P
end

function gradient(r::Array, alpha::Real, beta::Real, N::Integer)
	# function [dP] = basis.gradient(r, alpha, beta, N);
	# Purpose: Evaluate the derivative of the Jacobi polynomial of type (alpha,beta)>-1,
	# at points r for order N and returns dP[1:length(r))]
	r = reshape(r, size(r)[1], 1)
	dP = zeros(length(r), 1);
	if N == 0
		dP[:,:] = 0.0;
	else
		dP = sqrt(N * (N + alpha + beta + 1)) * polynomial(r[:], alpha + 1, beta + 1, N - 1);
	end
	return dP
end

function gauss(alpha::Real, beta::Real, N::Integer)
	# function [x,w] = basis.gauss(alpha,beta,N)
	# Purpose: Compute the N'th order Gauss quadrature points, x,
	# and weights, w, associated with the Jacobi
	# polynomial, of type (alpha,beta) > -1 ( <> -0.5).

	x = Float64[]
	w = Float64[]
	if N == 0
		push!(x, (alpha - beta) / (alpha + beta + 2))
		push!(w, 2)
		return x, w
	end

	# Form symmetric matrix from recurrence.
	# J = zeros(N+1);
	h1 = 2 * [0:N] + alpha + beta
	J = diagm(vec(-1 / 2 * (alpha^2 - beta^2) ./ (h1 + 2) ./ h1)) + diagm(2. / (h1[1:N] + 2) .* sqrt([1:N] .* ([1:N] + alpha + beta) .* ([1:N] + alpha) .* ([1:N] + beta) ./ (h1[1:N] + 1) ./ (h1[1:N] + 3)), 1)
	if alpha + beta < 10 * eps()
		J[1,1] = 0.0
	end
	J = J + J'
	# Compute quadrature by eigenvalue solve
    (D, V) = eig(J)
	x = reshape(D, size(D)[1], 1)
	w = (V[1,:]').^2 * 2^(alpha + beta + 1) / (alpha + beta + 1) * gamma(alpha + 1) * gamma(beta + 1) / gamma(alpha + beta + 1)
	return x, w
end

function gll(alpha::Real, beta::Real, N::Integer)
	# function [x] = basis.gll (alpha,beta,N)
	# Purpose: Compute the N'th order Gauss Lobatto quadrature
	# points, x, associated with the Jacobi polynomial,
	# of type (alpha,beta) > -1 ( <> -0.5).
	x = zeros(N + 1, 1)
	w = Float64[]
	if N == 1
		x[1] = -1.0
		x[2] = 1.0
		push!(w, 1.0)
		push!(w, 1.0)
		return x, w
	end
	xint, wq = gauss(alpha + 1, beta + 1, N - 2)
	x = hcat(-1, xint', 1)'
	# compute the weights
	w = polynomial(x, alpha, beta, N)
	adgammaN = (2.0 * N + alpha + beta + 1.0) / (N * (alpha + beta + N + 1.0))
	w = w .* w
	w = adgammaN ./ w
	w[1] = w[1] * (1.0 + alpha)
	w[end] = w[end] * (1.0 + beta)
	return x, w
end

end
