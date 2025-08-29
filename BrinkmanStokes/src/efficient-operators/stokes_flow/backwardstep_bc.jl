###BACKWARDSTEP_BC specifies streamfunction associated with backward step flow
# input
#   x       x boundary coordinate vector
#   y       y boundary coordinate vector
function backwardstep_bc(xbd, ybd)

    bc = 0 * xbd
    k = findall(xbd .== -1)
    bc[k] = 2ybd[k] .* ybd[k] .* (1 - 2ybd[k] / 3)
    k = findall(ybd .== 1)
    bc[k] = 2 / 3
    return bc

end
