using PyPlot

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

	# Create output directory for plots
	plots_dir = joinpath(@__DIR__, "../../../output/plots")
	mkpath(plots_dir)

	figure(figsize = (10, 8))
	streamplot(x, y, ux, uy, density = 4, color = ux)
	axis([-1, 1, -1, 1])
	title("Streamlines")
	colorbar()
	savefig(joinpath(plots_dir, "streamlines.png"), dpi = 150, bbox_inches = "tight")
	println("Streamlines plot saved to: $(joinpath(plots_dir, "streamlines.png"))")

	figure(figsize = (10, 8))
	# pcolor(x, y, kp, cmap = "Greys");
	quiver(x, y, ux, uy, ux, scale = 20)
	axis([-1, 1, -1, 1])
	title("Velocity Field")
	savefig(joinpath(plots_dir, "velocity_field.png"), dpi = 150, bbox_inches = "tight")
	println("Velocity field plot saved to: $(joinpath(plots_dir, "velocity_field.png"))")

	# Close figures to free memory
	close("all")
end
