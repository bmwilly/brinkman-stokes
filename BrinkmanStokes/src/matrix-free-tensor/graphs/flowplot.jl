using PyPlot
include("../../output_utils.jl")

function flowplot(sol, kparams, domain)

    x = kparams["x"]
    y = kparams["y"]
    xy = kparams["xy"]
    xyp = kparams["xyp"]
    mv = kparams["mv"]
    bound = kparams["bound"]
    ae = kparams["ae"]
    bxe = kparams["bxe"]
    bye = kparams["bye"]
    bbxe = kparams["bbxe"]
    bbye = kparams["bbye"]
    nvtx = length(xy[:, 1])
    nu = 2nvtx
    np = 3length(xyp[:, 1])

    u = sol[1:nu]
    p = sol[(nu + 1):end]
    ux = reshape(u[1:nvtx], length(x), length(y))'
    uy = reshape(u[(nvtx + 1):end], length(x), length(y))'

    # Get msize from kparams
    msize = kparams["msize"]

    # Create streamlines plot using unified output system
    streamlines_file = get_output_file("matrix-free-tensor", domain, msize, "streamlines.png"; subdir = "plots")
    figure(figsize = (10, 8))
    streamplot(x, y, ux, uy, density = 4, color = ux)
    axis([-1, 1, -1, 1])
    title("Streamlines (Domain=$domain, Size=$msize)")
    colorbar()
    savefig(streamlines_file, dpi = 150, bbox_inches = "tight")
    println("Streamlines plot saved to: $(streamlines_file)")

    # Create velocity field plot using unified output system
    velocity_file = get_output_file("matrix-free-tensor", domain, msize, "velocity_field.png"; subdir = "plots")
    figure(figsize = (10, 8))
    # pcolor(x, y, kp, cmap = "Greys");
    quiver(x, y, ux, uy, ux, scale = 20)
    axis([-1, 1, -1, 1])
    title("Velocity Field (Domain=$domain, Size=$msize)")
    savefig(velocity_file, dpi = 150, bbox_inches = "tight")
    println("Velocity field plot saved to: $(velocity_file)")

    # Close figures to free memory
    return close("all")
end
