include("regcavity_flow.jl")
include("poiseuille_flow.jl")
include("bfun.jl")

function gbc(domain, kparams)

	xy = kparams["xy"]
	xyp = kparams["xyp"]
	mv = kparams["mv"]
	bound = kparams["bound"]

	# get variables
	nvtx = length(xy[:, 1])
	nu = 2nvtx
	np = 3length(xyp[:, 1])
	nel = length(mv[:, 1])
	w = zeros(nu + np, 1)

	# get boundary
	xbd = xy[bound, 1]
	ybd = xy[bound, 2]
	if domain == 1
		(bcx, bcy) = regcavity_flow(xbd, ybd)
	elseif domain == 2
		(bcx, bcy) = poiseuille_flow(xbd, ybd)
	end

	# impose boundary conditions
	bccx = zeros(nvtx)
	bccx[bound] = bcx
	bccy = zeros(nvtx)
	bccy[bound] = bcy

	bc = [bccx; bccy; zeros(np, 1)]
	w -= bfun(bc, kparams)
	w = vec(w[1:np])
	return w
end
