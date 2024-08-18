###REGCAVITY_FLOW specifies regularized cavity flow boundary condition
# input
#   xbd         x coordinate vector
#   ybd         y coordinate vector
function regcavity_flow(xbd, ybd)

    bcx = 0 * xbd; bcy = 0 * xbd
    k = findall((ybd .== 1) .& (xbd .> -1) .& (xbd .< 1))
    bcx[k] = (1 .- xbd[k] .* xbd[k]) .* (1 .+ xbd[k] .* xbd[k])
    (bcx, bcy)

end
