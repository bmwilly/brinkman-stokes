using Dierckx
# using GridInterpolations
# using Gaston
using Gadfly

reload("stokes_flow/solve_stokes.jl")
reload("helpers/input.jl")
reload("stokes_flow/streambc.jl")
reload("grids/findobsXY.jl")

domain = int(input("Choose domain (1/square, 2/channel with obstacles): "))
# domain = 1
sol = solve_stokes(domain)
println("done")

xst = sol["xst"]; By = sol["By"]; Bx = sol["Bx"]; A = sol["A"];
xy = sol["xy"]; xyp = sol["xyp"]; x = sol["x"]; y = sol["y"];
bound = sol["bound"]; bndxy = sol["bndxy"]; bnde = sol["bnde"];
obs = sol["obs"]

nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
Asv = A[1:nvtx, 1:nvtx]

# compute auxiliary quantities
u = xst[1:nu]
p = xst[nu + 1:end]
f = [By -Bx] * u
(Asv, fsv) = streambc(Asv, f, xy, bound, domain)
phi = Asv\fsv
writecsv("/Users/bwilliams/Downloads/phij.csv", phi)

# grid = RectangleGrid(xy[:, 1], xy[:, 2])
# grid = RectangleGrid(vec(x), vec(y))
# interpolate(grid, phi, [x[1] y[1]])

# p = p[1:3:end]
# spl = Spline2D(
#   xyp[:, 1], xyp[:, 2], p;
#   kx = 4, ky = 5, s = length(xyp[:, 1])
# )
# zp = evalgrid(spl, vec(x), vec(y))

# plot(x = x, y = y, z = zp, Geom.contour)

# ## plot velocity
# xsol = unique(x)
# ysol = unique(y)
# # phi2 = vcat(zeros(length(xsol) * length(ysol) - length(phi), 1), phi)
# phi2 = vcat(phi, zeros(length(xsol) * length(ysol) - length(phi), 1))
# xysol = reshape(phi2, length(xsol), length(ysol))
#
# plot(x = xsol, y = ysol, z = xysol, Geom.contour)

# (X, Y) = meshgrid(x, y)
# spl = Spline2D(
#   xy[:, 1], xy[:, 2], phi;
#   kx = 4, ky = 4, s = length(xy[:, 1])
# )
# z = evalgrid(spl, vec(x), vec(y))
# # phi = vcat(zeros(length(x)))
# # phi = [phi[1:75]; zeros(9, 1); phi[76:end]]
# # phi1 = reshape(phi[1:63], 7, 9)'
# # phi2 = reshape(phi[64:64+36-1], 4, 9)'
# # phi3 = reshape(phi[64+36:end], 22, 9)'
# # z = [phi1 phi2 phi3]
# # z = reshape(xysol, length(x), length(y))
#
#
# # zlevels = -0.001:-0.001:-0.2
#
# if size(obs, 1) != 0
#   KK = findobsXY(obs, X, Y, bndxy)
#   z[KK] = NaN
# end
#
# writecsv("/Users/bwilliams/Downloads/zj.csv", z)
#
z = readcsv("/Users/bwilliams/Downloads/xysol.csv")
z = z'
writecsv("/Users/bwilliams/Downloads/z.csv",z)
plot(x = x, y = y, z = z, Geom.contour)


# plot(bndxy[bnde[1, 1], 1], bndxy[bnde[1, 2], 1], bndxy[bnde[1, 1], 2], bndxy[bnde[1, 2], 2])

# plot(x = [1 2 3], y = [4 5 6], Geom.point)
# plot(x = bndxy[:, 1], y = bndxy[:, 2], Geom.point)
#
# plot(
#   layer(x = bndxy[1:4, 1], y = bndxy[1:4, 2], Geom.path),
#   layer(x = bndxy[5:8, 1], y = bndxy[5:8, 2], Geom.path)
# )
# #
p = plot(
  layer(x = x, y = y, z = z, Geom.contour),
  layer(x = [bndxy[5:8, 1], bndxy[5, 1]], y = [bndxy[5:8, 2], bndxy[5, 2]], Geom.path)
)
draw(PNG("graphs/mat-flow.png", 8inch, 4inch), p)

# plotly attempt
# z = reshape(phi, length(x), length(y))
# data = [["x" => x, "y" => y, "z" => z, "type" => "contour"]]
# response = Plotly.plot(data)
# plot_url = response["url"]
