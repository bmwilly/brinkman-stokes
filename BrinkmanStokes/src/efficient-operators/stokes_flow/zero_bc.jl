###ZERO_BC specifies streamfunction associated with enclosed flow
# input
#   x       x boundary coordinate vector
#   y       y boundary coordinate vector
function zero_bc(xbd, ybd)
    return bc = 0 * xbd
end
