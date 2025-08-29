###HO_AFUN matrix-free A block operator for high order stokes
function ho_afun_nobc(u, params)
    mesh = params["mesh"]; order = params["order"]
    dof = params["dof"]; ne = params["ne"]; NP = params["NP"]
    bdy = params["bdy"]
    refel = params["refel"]
    centers = params["centers"]
    w = zeros(length(u))
    Ux = zeros(NP, ne); Uy = zeros(NP, ne)
    eMats = zeros(NP * ne, NP)

    # loop over elements
    for e in 1:ne
        idx = Mesh.get_node_indices(mesh, e, order)
        Ux[:, e] = u[idx]
        Uy[:, e] = u[idx + dof]
        pts = Mesh.element_nodes(mesh, e, refel)
        (detJac, Jac) = Mesh.geometric_factors(mesh, refel, pts)
        brinkman_pts = Mesh.brinkman_tensor(pts, centers)
        eMat = Mesh.element_stiffness_brinkman(mesh, e, refel, detJac, Jac, brinkman_pts)
        eMats[((e - 1) * NP + 1):(e * NP), :] = eMat
    end
    Wx = eMats * Ux; Wy = eMats * Uy
    for e in 1:ne
        idx = Mesh.get_node_indices(mesh, e, order)
        w[idx] += Wx[((e - 1) * NP + 1):(e * NP), e]
        w[idx + dof] += Wy[((e - 1) * NP + 1):(e * NP), e]
    end
    return vec(w)
end
