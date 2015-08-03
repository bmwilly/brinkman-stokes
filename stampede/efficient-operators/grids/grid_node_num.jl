# reload("grids/q2grid_node_num.jl")
# reload("grids/q3grid_node_num.jl")
# reload("grids/q4grid_node_num.jl")

###GRID_NODE_NUM gets number of nodes for grid made of specified element
# input
#   nelemx          number of elements along the x direction
#   nelemy          number of elements along the y direction
#   element         element type
# output
#   node_num        number of nodes
function grid_node_num(nelemx, nelemy, p)
  # if element == 2
  #   node_num = q2grid_node_num(nelemx, nelemy)
  # elseif element == 3
  #   node_num = q3grid_node_num(nelemx, nelemy)
  # elseif element == 4
  #   node_num = q4grid_node_num(nelemx, nelemy)
  # else
  #   error("invalide element type ...aborted...")
  # end
  # node_num
  (p*nelemx + 1) * (p*nelemy + 1)
end
