###Q2GRID_ELEMENT_NUM counts the elements in a Q2 element grid
# input
#   nelemx          number of elements along the x direction
#   nelemy          number of elements along the y direction
# output
#   element_num     number of elements in the grid
function q2grid_element_num(nelemx, nelemy)
  nelemx * nelemy
end
