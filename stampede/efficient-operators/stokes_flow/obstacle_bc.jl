###OBSTACLE_BC specifies streamfunction associated with obstacle flow
# input
#   xbd       x boundary coordinate vector
#   ybd       y boundary coordinate vector
function obstacle_bc(xbd, ybd)

  bc = 0 * xbd
  k = find(xbd .== 0); bc[k] = ybd[k] .^ 3/3
  k = find(ybd .== 1); bc[k] = 2/3
  k = find(ybd .== -1); bc[k] = -2/3
  bc
  
end
