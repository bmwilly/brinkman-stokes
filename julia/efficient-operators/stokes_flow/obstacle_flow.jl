###OBSTACLE_FLOW specifies flow over obstacle boundary condition condition
# input
#   xbd       x coordinate vector
#   ybd       y coordinate vector
function obstacle_flow(xbd, ybd)

  bcx = 0 * xbd; bcy = 0 * xbd;
  k = find(xbd .== 0)
  bcx[k] = (1 - ybd[k]) .* (1 + ybd[k])
  (bcx, bcy)

end
