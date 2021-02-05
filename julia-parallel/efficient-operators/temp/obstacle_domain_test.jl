reload("grids/obstacle_domain.jl")
obstacle_domain()

# println("Grid generation for domain with obstacle.")
# nc = input("grid parameter (3 for underlying 8x20 grid): ")
#
# if nc < 3
#   error("illegal parameter choice, try again")
# elseif nc == 3
#   println("Warning: increasing obstacle size for nc = 3")
# end
#
# n = 2^nc
# obs = []; xy = []; sbnde = [];
#
# # problem specific parameters
# # bndxy: vertices that define the domain and obstacle(s) (x-coord y-coord)
# # bnde: boundary edges (node1 node2)
# # obs: obstacles (node1 node2 node3 node4)
# # sbnde: boundary edges near which stretching is needed (edge1 edge2 ...)
# # "obs" and/or "sbnde" can be absent if there is no obstacle in the problem and/or only uniform grid is needed
#
# if nc == 3
#   bndxy = float([
#     0    -1
#     8     -1
#     8     1
#     0     1
#     1.5   -0.5
#     2.5   -0.5
#     2.5   0.5
#     1.5   0.5
#   ])
# else
#   bndxy = float([
#     0    -1
#     8     -1
#     8     1
#     0     1
#     1.75  -0.25
#     2.25  -0.25
#     2.25  0.25
#     1.75  0.25
#   ])
# end
#
# bnde = [
#   1 2 1
#   2 3 0
#   3 4 1
#   4 1 1
#   5 6 1
#   6 7 1
#   7 8 1
#   8 5 1
# ]
#
# obs = [5 6 7 8]
#
# sbnde = [5 6 7 8]
#
# # compute mesh size h and uniform x, y coordinate (problem nonspecific)
# h = min(maximum(bndxy[:, 2]) - minimum(bndxy[:, 2]), maximum(bndxy[:, 1]) - minimum(bndxy[:, 1])) / n
# x = [minimum(bndxy[:, 1]):h:maximum(bndxy[:, 1])]
# y = [minimum(bndxy[:, 2]):h:maximum(bndxy[:, 2])]
#
# # generate xy coordinates and element mapping matrix (problem specific)
# (xy1, mv1) = rec_domain_q2(0, bndxy[obs[1], 1], -1, 1, h, xy)
# xy = xy1
# mv = mv1
#
# # (xy2, mv2) = rec_domain_q2(bndxy[obs[1], 1], bndxy[obs[2], 1], -1, bndxy[obs[1], 2], h, xy)
# x1 = bndxy[obs[1],1]; x2 = bndxy[obs[2],1]; y1 = -1; y2 = bndxy[obs[1],2]; xyo = xy
#
# xmin = min(x1, x2)
# xmax = max(x1, x2)
# ymin = min(y1, y2)
# ymax = max(y1, y2)
#
# x = [xmin:h:xmax]
# y = [ymin:h:ymax]
# nx = length(x) - 1
# ny = length(y) - 1
# nxyo = size(xyo, 1)
# nvtx = (nx + 1) * (ny + 1)
# (X, Y) = meshgrid(x, y)
# xx = reshape(X', nvtx, 1)
# yy = reshape(Y', nvtx, 1)
# xy = [xx[:] yy[:]]
# xyg = [1:nvtx] + nxyo
#
# if nxyo != 0
#   for i = 1:nvtx
#     ix = findall((xyo[:, 1] .== xy[i, 1]) & (xyo[:, 2] .== xy[i, 2]))
#     if size(ix, 1) != 0
#       xyg[i] = ix[1]
#       if i < nvtx
#         xyg[(i+1):nvtx] -= 1
#       end
#     end
#     if xyg[i] > nxyo
#       if xyn == []
#         xyn = xy[i, :]
#       else
#         xyn = [xyn; xy[i, :]]
#       end
#     end
#   end
#   xy = xyn
# end
#
# kx = 1; ky = 1; mel = 0;
# nvv = zeros(1, 9);
# mv = zeros(int(nx/2 * ny/2), 9)
# for j = 1:ny/2
#   for i = 1:nx/2
#     mref = (nx + 1) * (ky - 1) + kx
#     mel +=1
#     nvv[1] = xyg[mref]
#     nvv[2] = xyg[mref + 2]
#     nvv[3] = xyg[mref + 2nx + 4]
#     nvv[4] = xyg[mref + 2nx + 2]
#     nvv[5] = xyg[mref + 1]
#     nvv[6] = xyg[mref + nx + 3]
#     nvv[7] = xyg[mref + 2nx + 3]
#     nvv[8] = xyg[mref + nx + 1]
#     nvv[9] = xyg[mref + nx + 2]
#     mv[mel; 1:9] = nvv[1:9]
#     kx += 2
#   end
#   ky += 2
#   kx = 1
# end
#
# (xy, mv)
