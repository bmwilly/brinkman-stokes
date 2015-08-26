reload("helpers/input.jl")
reload("helpers/meshgrid.jl")

###CAVITY_DOMAIN square cavity Q2 grid generator
function cavity_domain(msize)

  ## define geometry
  println("Grid generation for cavity domain.")
  n = 2^msize
  np = int(n/2)
  nel = int(np^2)

  # y-direction
  yy = [1/np:1/np:1]
  ypos = [0, yy]
  yneg = -yy[length(yy):-1:1]
  y = [yneg, ypos]'
  x = y

  # compute biquadratic element coordinates
  nvtx = (n+1) * (n+1)
  (X,Y) = meshgrid(x,y)
  xx = reshape(X', nvtx, 1)
  yy = reshape(Y', nvtx, 1)
  xy = [xx[:] yy[:]]

  kx = 1
  ky = 1
  mel = 0
  mv = zeros(Int64, nel, 9)
  for j = 1:np
      for i = 1:np
          mref = (n+1)*(ky-1) + kx
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
          mv[mel, 1:9] = nvv[1:9]
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

  ## specify boundary information for graphics
  # bndxy: (x,y)-coordinates of vertices that define the domain and obstacle(s)
  # bnde: boundary edges (node 1 node 2)
  # obs: obstacles (node1 node2 node3 node4)
  # sbnde: boundary edges near which stretching is needed (edge1 edge2 ...)
  # "obs" and/or "sbnde" can be absent if there is no obstacle in the problem
  bndxy = [-1 -1; 1 -1; 1 1; -1 1]
  bnde = [1 2 1; 2 3 1; 3 4 1; 4 1 1]
  obs = []
  sbnde = [1 2 3 4]

  cavity_grid = {
    "mv" => mv,
    "xy" => xy,
    "bound" => bound,
    "mbound" => mbound,
    "x" => x,
    "y" => y,
    "bndxy" => bndxy,
    "bnde" => bnde,
    "obs" => obs
  }
end
