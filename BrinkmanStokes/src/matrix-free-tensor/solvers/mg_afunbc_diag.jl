###MG_AFUNBC_DIAG modified diagonal operator with zero boundary condition imposed
function mg_afunbc_diag(u, mparams)

    ev = mparams["ev"]
    bound = mparams["bound"]
    ae = mparams["ae"]
    nel = length(ev[:, 1])
    aes = squeeze(ae[1, :, :], 1)
    w = zeros(length(u), 1)

    for e in 1:nel
        ind = ev[e, :]'
        indbd = findall(in(bound), ind)
        indint = setdiff(int(linspace(1, 4, 4)), indbd)
        indb = ind[indbd]

        ue = u[ind]
        ue[indbd] = zeros(1, length(indbd))

        we = diagm(diag(aes)) * ue

        # impose boundary condition
        we[indbd] = u[indb]

        # add interior points
        w[ind[indint]] += we[indint]

        # add boundary points only when entry == 0
        wb = w[indb]
        ind0 = findall(wb .== 0)
        web = we[indbd]
        w[indb[ind0]] += web[ind0]

    end

    return vec(w)

end
