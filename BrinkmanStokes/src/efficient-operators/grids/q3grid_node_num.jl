###Q3GRID_NODE_NUM counts the nodes in a Q3 element grid
# input
#   nelemx          number of elements along the x direction
#   nelemy          number of elements along the y direction
# output
#   node_num        number of nodes in the grid
function q3grid_node_num(nelemx, nelemy)
    (3nelemx + 1) * (3nelemy + 1)
end
