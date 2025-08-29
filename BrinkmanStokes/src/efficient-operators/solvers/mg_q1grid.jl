###MG_Q1GRID bilinear element grid generator for GMG
# input
# 	x 				x-grid coordinates
# 	y 				y-grid coordinates
# 	xy 				vertex coordinate vector
# 	mv 				Q2 macroelement mapping matrix
# 	bound 			boundary vertex vector
# 	mbound 			macroelement boundary edge vector
# output
# 	ev 				element mapping matrix
# 	ebound 			element boundary mapping matrix
function mg_q1grid(x, y, xy, mv, bound, mbound)

	xx = xy[:, 1]
	yy = xy[:, 2]
	nvtx = length(xx)
	adj = spzeros(nvtx, nvtx)
	mel = length(mv[:, 1])
	nel = 4mel
	ev = zeros(Int, nel, 4)
	# Calculate total boundary elements (each boundary element creates 2 entries)
	total_boundary_elements = 2 * (length(findall(mbound[:, 2] .== 1)) +
								   length(findall(mbound[:, 2] .== 2)) +
								   length(findall(mbound[:, 2] .== 3)) +
								   length(findall(mbound[:, 2] .== 4)))
	ebound = zeros(total_boundary_elements, 2)

	# loop over macroelements
	for k in 1:mel

		# first element
		ke = 4k - 3
		ev[ke, 1] = Int(mv[k, 1])
		ev[ke, 2] = Int(mv[k, 5])
		ev[ke, 3] = Int(mv[k, 9])
		ev[ke, 4] = Int(mv[k, 8])

		# second element
		ke = 4k - 2
		ev[ke, 1] = Int(mv[k, 5])
		ev[ke, 2] = Int(mv[k, 2])
		ev[ke, 3] = Int(mv[k, 6])
		ev[ke, 4] = Int(mv[k, 9])

		# third element
		ke = 4k - 1
		ev[ke, 1] = Int(mv[k, 9])
		ev[ke, 2] = Int(mv[k, 6])
		ev[ke, 3] = Int(mv[k, 3])
		ev[ke, 4] = Int(mv[k, 7])

		# fourth element
		ke = 4k
		ev[ke, 1] = Int(mv[k, 8])
		ev[ke, 2] = Int(mv[k, 9])
		ev[ke, 3] = Int(mv[k, 7])
		ev[ke, 4] = Int(mv[k, 4])

	end

	# define element edges
	ect = 1
	# bottom boundary edges
	k1 = findall(mbound[:, 2] .== 1)
	for idx in k1
		k                = mbound[idx, 1]
		ebound[ect, 1]   = 4k - 3
		ebound[ect+1, 1] = 4k - 2
		ebound[ect, 2]   = 1
		ebound[ect+1, 2] = 1
		ect              += 2
	end

	# right boundary edges
	k2 = findall(mbound[:, 2] .== 2)
	for idx in k2
		k                = mbound[idx, 1]
		ebound[ect, 1]   = 4k - 2
		ebound[ect+1, 1] = 4k - 1
		ebound[ect, 2]   = 2
		ebound[ect+1, 2] = 2
		ect              += 2
	end

	# top boundary edges
	k3 = findall(mbound[:, 2] .== 3)
	for idx in k3
		k                = mbound[idx, 1]
		ebound[ect, 1]   = 4k - 1
		ebound[ect+1, 1] = 4k
		ebound[ect, 2]   = 3
		ebound[ect+1, 2] = 3
		ect              += 2
	end

	# left boundary edges
	k4 = findall(mbound[:, 2] .== 4)
	for idx in k4
		k                = mbound[idx, 1]
		ebound[ect, 1]   = 4k
		ebound[ect+1, 1] = 4k - 3
		ebound[ect, 2]   = 4
		ebound[ect+1, 2] = 4
		ect              += 2
	end

	return (ev, ebound)

end
