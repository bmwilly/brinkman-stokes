reload("diffusion/hshape.jl")

###HDERIV evaluates derivatives of arbitrary order shape functions
#input
#    s           reference element x coordinate
#    t           reference element y coordinate
#    xl          physical element x vertex coordinates
#    yl          physical element y vertex coordinates
#output
#    psi         element-wise shape functions
#    dpsidx      x derivatives of psi
#    dpsidy      y derivatives of psi
function hderiv(s, t, xl, yl)

    nel = length(xl[:, 1])
    ninel = length(xl[1, :])
    p = Int(sqrt(ninel))
    nn = (p + 1) * (p + 1)
    zerov = zeros(nel, 1)
    onev = ones(nel, 1)

    # evaluate shape functions
    (psie, dpsids, dpsidt) = hshape(s, t, p)

    # local derivatives
    dxds = copy(zerov)
    dxdt = copy(zerov)
    dyds = copy(zerov)
    dydt = copy(zerov)
    jac = copy(zerov)
    invjac = copy(zerov)

    for ivtx in 1:ninel
        dxds[:] += xl[:, ivtx] .* onev * dpsids[ivtx]
        dxdt[:] += xl[:, ivtx] .* onev * dpsidt[ivtx]
        dyds[:] += yl[:, ivtx] .* onev * dpsids[ivtx]
        dydt[:] += yl[:, ivtx] .* onev * dpsidt[ivtx]
    end

    jac[:] = dxds[:] .* dydt[:] - dxdt[:] .* dyds[:]

    # check element Jacobian
    # if any(jac .< 1e-9)
    #     println("Bad element warning ...")
    #     if any(jac .<= 0)
    #         error("singular Jacobian ... aborted ...")
    #     end
    # end

    invjac[:] = onev ./ jac[:]

    psi = zeros(nel, nn)
    dpsidx = zeros(nel, nn)
    dpsidy = zeros(nel, nn)

    for ivtx in 1:nn
        psi[:, ivtx] = psie[ivtx] * onev
        dpsidx[:, ivtx] = dpsids[ivtx] .* dydt[:] - dpsidt[ivtx] .* dyds[:]
        dpsidy[:, ivtx] = -dpsids[ivtx] .* dxdt[:] + dpsidt[ivtx] .* dxds[:]
    end

    return (jac, invjac, psi, dpsidx, dpsidy)

end
