###GRID_ELEMENT_NUM counts the elements in a Q_ element grid
# input
#   nelemx          number of elements along the x direction
#   nelemy          number of elements along the y direction
# output
#   element_num     number of elements in the grid
function grid_element_num(nelemx, nelemy)
    nelemx * nelemy
end
