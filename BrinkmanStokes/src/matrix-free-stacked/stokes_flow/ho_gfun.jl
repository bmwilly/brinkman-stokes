###HO_GFUN matrix-free mass matrix
function ho_gfun(u, params)
    mesh = params["mesh"]; order = params["order"]
    dof = params["dof"]
    bdy = params["bdy"]
    refel = params["refel"]
    w = zeros(length(u))

    # loop over elements
    for e in 1:ne
        idx = Mesh.get_node_indices(mesh, e, order)
        pts = Mesh.element_nodes(mesh, e, refel)
        (detJac, Jac) = Mesh.geometric_factors(mesh, refel, pts)
        brinkman_pts = Mesh.brinkman_tensor(pts, centers)
        eMat = Mesh.element_mass_brinkman(mesh, e, refel, detJac, Jac, brinkman_pts)
        w[idx] += eMat * u[idx]
        w[idx + dof] += eMat * u[idx + dof]
    end
    return vec(w)
end
