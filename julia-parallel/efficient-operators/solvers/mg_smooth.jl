###MG_SMOOTH smoothers for GMG on a square domain
# input
# 	As 				array of dicts containing coefficient matrices
# 	level 			finest level of grid
# 	sweeps			number of directions used for Gauss-Seidel
# 	smooth 			choice of smoother (1:Jacobi, 2:Gauss-Seidel, 3:ILU)
# 	stype 			type of smoother (1:point, 2:line)
# output
# 	Qs 				dict containing smoothing operator in factored form
function mg_smooth(As, level, sweeps, smooth, stype)

	Qs = Array(Dict, Int(level))

	for i = level:-1:2
		A = As[i]["matrix"]
		N = size(A, 1)
		n = sqrt(N)

		# line Gauss-Seidel
		if stype == 2
			Q1 = tril(A, 1)
			(L1, U1) = lu(Q1)

			if sweeps >= 2
				Q2 = diagm(diag(A, 0)) +
					diagm(diag(A, -n), -n) +
					diagm(diag(A, n), n) +
					diagm(diag(A, -n-1), -n-1) +
					diagm(diag(A, -1), -1) +
					diagm(diag(A, n-1), n-1)

				(L2, U2) = lu(Q2)

				if sweeps >= 3
					Q3 = triu(A, -1)
					(L3, U3) = lu(Q3)

					if sweeps == 4
						Q4 = diagm(diag(A, 0)) +
							diagm(diag(A, -n), -n) +
							diagm(diag(A, n), n) +
							diagm(diag(A, n+1), n+1) +
							diagm(diag(A, 1), 1) +
							diagm(diag(A, -n+1), -n+1)

						(L4, U4) = lu(Q4)

					else
						L4 = spzeros(N, N)
						U4 = L4
					end
				else
					L3 = spzeros(N, N)
					U3 = L3; L4 = L3; U4 = L3
				end
			else
				L2 = sparse(N, N)
				U2 = L2; L3 = L2; U3 = L2; L4 = L2; U4 = L2
			end

		# point smoothers
		else

			#ILU
			if smooth == 3
				#TODO

			# point Gauss-Seidel
			elseif smooth == 2
				Q1 = tril(A, 0)
				(L1, Q1) = lu(Q1)
				l1, l2 = size(L1)
				L2 = spzeros(l1, l2);
				U2 = L2; L3 = L2; U3 = L2; L4 = L2; U4 = L2;

			# point damped Jacobi
			else
				omega = 8/9 # relaxation factor
				Q1 = (1/omega) * spdiagm(diag(A), 0, N, N)
				(L1, U1) = lu(full(Q1))
				L1 = sparse(L1); U1 = sparse(U1)
				l1, l2 = size(L1)
				L2 = spzeros(l1, l2)
				U2 = L2; L3 = L2; U3 = L2; L4 = L2; U4 = L2
			end
		end

		qs = Dict(
			"L1" => L1, "L2" => L2, "L3" => L3, "L4" => L4,
			"U1" => U1, "U2" => U2, "U3" => U3, "U4" => U4
		)
		Qs[i] = qs

	end # for loop

	Qs

end
