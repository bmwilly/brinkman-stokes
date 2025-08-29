reload("diffusion/q2shape.jl")
reload("diffusion/q3shape.jl")

###Q3DERIV evaluates derivatives of bicubic shape functions
#input
#    s           reference element x coordinate
#    t           reference element y coordinate
#    xl          physical element x vertex coordinates
#    yl          physical element y vertex coordinates
#output
#    chi         element-wise shape functions
#    dchidx      x derivatives of chi
#    dchidy      y derivatives of chi
function q3deriv(s, t, xl, yl)

    nel = length(xl[:, 1])
    # zerov = zeros(length(xl[1, :]), 1)
    zerov = zeros(nel, 1)
    onev = ones(nel, 1)

    # evaluate bilinear shape functions
    (phie, dphids, dphidt) = shape(s, t)

    # evaluate biquadratic shape functions
    (psie, dpsids, dpsidt) = q2shape(s, t)

    # evaluate bicubic shape functions
    (chie, dchids, dchidt) = q3shape(s, t)

    # local derivatives
    dxds = copy(zerov)
    dxdt = copy(zerov)
    dyds = copy(zerov)
    dydt = copy(zerov)
    jac = copy(zerov)
    invjac = copy(zerov)

    # for ivtx = 1:4
    #     dxds[:] += xl[:,ivtx] .* onev * dphids[ivtx]
    #     dxdt[:] += xl[:,ivtx] .* onev * dphidt[ivtx]
    #     dyds[:] += yl[:,ivtx] .* onev * dphids[ivtx]
    #     dydt[:] += yl[:,ivtx] .* onev * dphidt[ivtx]
    # end

    for ivtx in 1:9
        dxds[:] += xl[:, ivtx] .* onev * dchids[ivtx]
        dxdt[:] += xl[:, ivtx] .* onev * dchidt[ivtx]
        dyds[:] += yl[:, ivtx] .* onev * dchids[ivtx]
        dydt[:] += yl[:, ivtx] .* onev * dchidt[ivtx]
    end

    # for ivtx = 1:16
    #     dxds[:] += xl[:,ivtx] .* onev * dchids[ivtx]
    #     dxdt[:] += xl[:,ivtx] .* onev * dchidt[ivtx]
    #     dyds[:] += yl[:,ivtx] .* onev * dchids[ivtx]
    #     dydt[:] += yl[:,ivtx] .* onev * dchidt[ivtx]
    # end

    jac[:] = dxds[:] .* dydt[:] - dxdt[:] .* dyds[:]

    # check element Jacobian
    # if any(jac .< 1e-9)
    #     println("Bad element warning ...")
    #     if any(jac .<= 0)
    #         error("singular Jacobian ... aborted ...")
    #     end
    # end

    invjac[:] = onev ./ jac[:]

    chi = zeros(nel, 16)
    dchidx = zeros(nel, 16)
    dchidy = zeros(nel, 16)

    for ivtx in 1:16
        chi[:, ivtx] = chie[ivtx] * onev
        dchidx[:, ivtx] = dchids[ivtx] .* dydt[:] - dchidt[ivtx] .* dyds[:]
        dchidy[:, ivtx] = -dchids[ivtx] .* dxdt[:] + dchidt[ivtx] .* dxds[:]
    end

    return (jac, invjac, chi, dchidx, dchidy)

end
