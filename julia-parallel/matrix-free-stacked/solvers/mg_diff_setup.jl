reload("solvers/mg_q1grid.jl")
reload("solvers/mg_q1diff.jl")
# reload("solvers/mg_zerobc.jl")
reload("solvers/mg_afunbc.jl")

###MG_DIFF_SETUP GMG diffusion problem on square domain
# input
# 	x 				x coordinate vector for coarse grid
# 	y 				y coordinate vector for coarse grid
# output
# 	Agal 			discrete diffusion operator
function mg_diff_setup(x, y)

	n = length(x) - 1; np = n/2; nq = n/4
  nel = int(np^2); nvtx = (n + 1) * (n + 1)
	(X, Y) = meshgrid(x, y)
	xx = reshape(X', nvtx, 1)
	yy = reshape(Y', nvtx, 1)
	xy = [xx[:] yy[:]]

	# assemply process
	kx = 1; ky = 1; mel = 0
	mv = zeros(Int64, nel, 9)
    mp = zeros(Int64, nel, 4)
    for j = 1:np
        for i = 1:np
            mref = (n + 1) * (ky - 1) + kx
            pref = (np + 1) * (j - 1) + i
            mel += 1
            nvv = zeros(9)
            nvv[1] = mref
            nvv[2] = mref + 2
            nvv[3] = mref + 2n + 4
            nvv[4] = mref + 2n + 2
            nvv[5] = mref + 1
            nvv[6] = mref + n + 3
            nvv[7] = mref + 2n + 3
            nvv[8] = mref + n + 1
            nvv[9] = mref + n + 2
            npp = zeros(4)
            npp[1] = pref
            npp[2] = pref + 1
            npp[3] = pref + np + 2
            npp[4] = pref + np + 1
            mv[mel, 1:9] = nvv[1:9]
            mp[mel, 1:4] = npp[1:4]
            kx += 2
        end
        ky += 2
        kx = 1
    end

    # compute boundary vertices and edges
    # four boundary edges
    k1 = find(xy[:,2] .== -1)
    e1 = Int[]
    for k = 1:mel
        if any(mv[k,5] .== k1)
						push!(e1, k)
        end
    end
    ef1 = ones(size(e1))

    k2 = find((xy[:,1] .== 1) & (xy[:,2] .< 1) & (xy[:,2] .> -1))
    e2 = Int[]
    for k = 1:mel
        if any(mv[k,6] .== k2)
						push!(e2, k)
        end
    end
    ef2 = 2*ones(size(e2))

    k3 = find(xy[:,2] .== 1)
    e3 = Int[]
    for k = 1:mel
        if any(mv[k,7] .== k3)
						push!(e3, k)
        end
    end
    ef3 = 3*ones(size(e3))

    k4 = find((xy[:,1] .== -1) & (xy[:,2] .< 1) & (xy[:,2] .> -1))
    e4 = Int[]
    for k = 1:mel
        if any(mv[k,8] .== k4)
						push!(e4, k)
        end
    end
    ef4 = 4*ones(size(e4))

    bound = sort([k1; k2; k3; k4])
    mbound = [e1' ef1'; e2' ef2'; e3' ef3'; e4' ef4']

    # set up matrices for Q1 approximation
    (ev, ebound) = mg_q1grid(x, y, xy, mv, bound, mbound)
    ae, f = mg_q1diff(xy, ev)
		mparams = {"xy" => xy, "ev" => ev, "bound" => bound, "ae" => ae}

end
