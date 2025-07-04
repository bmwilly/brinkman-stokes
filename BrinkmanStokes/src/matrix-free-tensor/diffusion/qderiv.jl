include("qshape.jl")

###QDERIV evaluates derivatives of biquadratic shape functions
#input
#    s           reference element x coordinate
#    t           reference element y coordinate
#    xl          physical element x vertex coordinates
#    yl          physical element y vertex coordinates
#output
#    psi         element-wise shape functions
#    dpsidx      x derivatives of psi
#    dpsidy      y derivatives of psi
function qderiv(s, t, xl, yl)

    nel = length(xl[:, 1])
    zerov = zeros(nel, 1)
    onev = ones(nel, 1)

    # evaluate bilinear shape functions
    (phie, dphids, dphidt) = shape(s, t)

    # evaluate biquadratic shape functions
    (psie, dpsids, dpsidt) = qshape(s, t)

    # local derivatives
    dxds = copy(zerov)
    dxdt = copy(zerov)
    dyds = copy(zerov)
    dydt = copy(zerov)
    jac = copy(zerov)
    invjac = copy(zerov)

    for ivtx = 1:4
        dxds[:] += xl[:,ivtx] .* onev * dphids[ivtx]
        dxdt[:] += xl[:,ivtx] .* onev * dphidt[ivtx]
        dyds[:] += yl[:,ivtx] .* onev * dphids[ivtx]
        dydt[:] += yl[:,ivtx] .* onev * dphidt[ivtx]
    end

    jac[:] = dxds[:] .* dydt[:] - dxdt[:] .* dyds[:]

    # check element Jacobian
    if any(jac .< 1e-9)
        println("Bad element warning ...")
        if any(jac .<= 0)
            error("singular Jacobian ... aborted ...")
        end
    end

    invjac[:] = onev ./ jac[:]

    psi = zeros(nel, 9)
    dpsidx = zeros(nel, 9)
    dpsidy = zeros(nel, 9)

    for ivtx = 1:9
        psi[:, ivtx] = psie[ivtx] * onev
        dpsidx[:, ivtx] = dpsids[ivtx] .* dydt[:] - dpsidt[ivtx] .* dyds[:]
        dpsidy[:, ivtx] = -dpsids[ivtx] .* dxdt[:] + dpsidt[ivtx] .* dxds[:]
    end

    (psi, dpsidx, dpsidy)

end
