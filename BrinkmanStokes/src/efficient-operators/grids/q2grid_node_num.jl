###Q2GRID_NODE_NUM counts the nodes in a Q2 element grid
# input
#   nelemx          number of elements along the x direction
#   nelemy          number of elements along the y direction
# output
#   node_num        number of nodes in the grid
function q2grid_node_num(nelemx, nelemy)
    (2nelemx + 1) * (2nelemy + 1)
end
