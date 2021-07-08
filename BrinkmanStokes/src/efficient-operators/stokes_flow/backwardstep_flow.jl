###BACKWARDSTEP_FLOW specifies backward step flow boundary condition
# input
#   xbd       x coordinate vector
#   ybd       y coordinate vector
function backwardstep_flow(xbd, ybd)

    bcx = 0 * xbd; bcy = 0 * xbd;
    k = findall(xbd .== -1)
    bcx[k] = 4ybd[k] .* (1 - ybd[k])
    (bcx, bcy)

end
