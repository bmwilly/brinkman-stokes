using PyPlot
include("../../output_utils.jl")

function flowplot(sol, domain)

	xst = sol["xst"]
	By = sol["By"]
	Bx = sol["Bx"]
	A = sol["A"]
	xy = sol["xy"]
	xyp = sol["xyp"]
	x = Array(sol["x"])
	y = Array(sol["y"])
	bound = sol["bound"]
	bndxy = sol["bndxy"]
	bnde = sol["bnde"]
	obs = sol["obs"]
	if domain == 3
		kappa = sol["kappa"]
		kp = reshape(kappa, length(x), length(y))'
	end

	nvtx = length(xy[:, 1])
	nu = 2nvtx
	np = 3length(xyp[:, 1])
	Asv = A[1:nvtx, 1:nvtx]
	x = vec(x)
	y = vec(y)
	xp = unique(xyp[:, 1])
	yp = unique(xyp[:, 2])

	# compute auxilliary quantities
	u = xst[1:nu]
	p = xst[nu+1:end]

	## plot velocity
	ux = reshape(u[1:nvtx], length(x), length(y))'
	uy = reshape(u[nvtx+1:end], length(x), length(y))'

	# Get msize from solution data
	msize = sol["msize"]

	# Create streamlines plot using unified output system
	streamlines_file = get_output_file("efficient-operators", domain, msize, "streamlines.png"; subdir = "plots")
	figure(figsize = (10, 8))
	streamplot(x, y, ux, uy, density = 4, color = ux)
	axis([-1, 1, -1, 1])
	title("Streamlines (Domain=$domain, Size=$msize)")
	colorbar()
	savefig(streamlines_file, dpi = 150, bbox_inches = "tight")
	println("Streamlines plot saved to: $(streamlines_file)")

	# Create velocity field plot using unified output system
	velocity_file = get_output_file("efficient-operators", domain, msize, "velocity_field.png"; subdir = "plots")
	figure(figsize = (10, 8))
	# pcolor(x, y, kp, cmap = "Greys");
	quiver(x, y, ux, uy, ux, scale = 20)
	axis([-1, 1, -1, 1])
	title("Velocity Field (Domain=$domain, Size=$msize)")
	savefig(velocity_file, dpi = 150, bbox_inches = "tight")
	println("Velocity field plot saved to: $(velocity_file)")

	# Close figures to free memory
	close("all")
end
