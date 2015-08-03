# include("../helpers/input.jl")
reload("helpers/input.jl")
reload("grids/rec_domain_q2.jl")
reload("grids/findboundary.jl")

###OBSTACLE_DOMAIN obstacle domain Q2 grid generator
function obstacle_domain()

  println("Grid generation for domain with obstacle.")
  nc = int(input("grid parameter (3 for underlying 8x20 grid): "))
  # nc = 3;

  if nc < 3
    error("illegal parameter choice, try again")
  elseif nc == 3
    println("Warning: increasing obstacle size for nc = 3")
  end

  n = 2^nc
  obs = []; xy = []; sbnde = [];

  # problem specific parameters
  # bndxy: vertices that define the domain and obstacle(s) (x-coord y-coord)
  # bnde: boundary edges (node1 node2)
  # obs: obstacles (node1 node2 node3 node4)
  # sbnde: boundary edges near which stretching is needed (edge1 edge2 ...)
  # "obs" and/or "sbnde" can be absent if there is no obstacle in the problem and/or only uniform grid is needed

  if nc == 3
    bndxy = float([
      0    -1
      8     -1
      8     1
      0     1
      1.5   -0.5
      2.5   -0.5
      2.5   0.5
      1.5   0.5
    ])
  else
    bndxy = float([
      0    -1
      8     -1
      8     1
      0     1
      1.75  -0.25
      2.25  -0.25
      2.25  0.25
      1.75  0.25
    ])
  end

  bnde = [
    1 2 1
    2 3 0
    3 4 1
    4 1 1
    5 6 1
    6 7 1
    7 8 1
    8 5 1
  ]

  obs = [5 6 7 8]

  sbnde = [5 6 7 8]

  # compute mesh size h and uniform x, y coordinate (problem nonspecific)
  h = min(maximum(bndxy[:, 2]) - minimum(bndxy[:, 2]), maximum(bndxy[:, 1]) - minimum(bndxy[:, 1])) / n
  x = [minimum(bndxy[:, 1]):h:maximum(bndxy[:, 1])]
  y = [minimum(bndxy[:, 2]):h:maximum(bndxy[:, 2])]

  # generate xy coordinates and element mapping matrix (problem specific)
  (xy1, mv1) = rec_domain_q2(0, bndxy[obs[1], 1], -1, 1, h, xy)
  xy = xy1
  mv = mv1

  (xy2, mv2) = rec_domain_q2(bndxy[obs[1], 1], bndxy[obs[2], 1], -1, bndxy[obs[1], 2], h, xy)
  xy = [xy; xy2]
  mv = [mv; mv2]

  (xy3, mv3) = rec_domain_q2(bndxy[obs[1], 1], bndxy[obs[2], 1], bndxy[obs[3], 2], 1, h, xy)
  xy = [xy; xy3]
  mv = [mv; mv3]

  (xy4, mv4) = rec_domain_q2(bndxy[obs[2], 1], 8, -1, 1, h, xy)
  xy = [xy; xy4]
  mv = [mv; mv4]

  # compute boundary vertices and edges (problem nonspecific)
  (bound, mbound) = findboundary(bndxy, bnde, xy, mv)

  obstacle_grid = {
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
