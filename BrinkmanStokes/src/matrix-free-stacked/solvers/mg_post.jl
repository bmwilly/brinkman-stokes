###MG_POST postsmoothing for GMG
# input
# 	A 				coefficient matrix
# 	xs0 			initial iterate
# 	f 				rhs
# 	ns 				number of smoothing steps
# 	Qs 				structure containing smoothing operators
# 	level 			grid level
# 	sweeps 			type of sweeping strategy used for Gauss-Seidel smoothing
# output
# 	xs 				result of postsmoothing applied to x0
function mg_post(A, params, xs0, f, ns, Qs, level, sweeps)

	xs = xs0
	if ns == 0
		return xs
	end

	L1 = Qs[level]["L1"];	U1 = Qs[level]["U1"]
	L2 = Qs[level]["L2"];	U2 = Qs[level]["U2"]
	L3 = Qs[level]["L3"];	U3 = Qs[level]["U3"]
	L4 = Qs[level]["L4"];	U4 = Qs[level]["U4"]

	r = f - A(xs, params)
	if sweeps >= 4
		r1,flag = gmres(u -> L4(u, params), r, length(r))
		r2,flag = gmres(u -> U4(u, params), r1, length(r1))
		xs += r2
		r = f - A(xs, params)
	end
	if sweeps >= 3
		r1,flag = gmres(u -> L3(u, params), r, length(r))
		r2,flag = gmres(u -> U3(u, params), r1, length(r1))
		xs += r2
		r = f - A(xs, params)
	end
	if sweeps >= 2
		r1,flag = gmres(u -> L2(u, params), r, length(r))
		r2,flag = gmres(u -> U2(u, params), r1, length(r1))
		xs += r2
		r = f - A(xs, params)
	end
	if sweeps >= 1
		r1,flag = gmres(u -> L1(u, params), r, length(r))
		r2,flag = gmres(u -> U1(u, params), r1, length(r1))
		xs += r2
	end

	if ns > 1
		k = 1
		while k < ns
			r = f - A(xs, params)
			if sweeps >= 4
				r1,flag = gmres(u -> L1(u, params), r, length(r))
				r2,flag = gmres(u -> U1(u, params), r1, length(r1))
				xs += r2
				r = f - A(xs, params)
			end
			if sweeps >= 3
				r1,flag = gmres(u -> L1(u, params), r, length(r))
				r2,flag = gmres(u -> U1(u, params), r1, length(r1))
				xs += r2
				r = f - A(xs, params)
			end
			if sweeps >= 2
				r1,flag = gmres(u -> L1(u, params), r, length(r))
				r2,flag = gmres(u -> U1(u, params), r1, length(r1))
				xs += r2
				r = f - A(xs, params)
			end
			if sweeps >= 1
				r1,flag = gmres(u -> L1(u, params), r, length(r))
				r2,flag = gmres(u -> U1(u, params), r1, length(r1))
				xs += r2
			end
			k += 1
		end
	end

	xs
end
