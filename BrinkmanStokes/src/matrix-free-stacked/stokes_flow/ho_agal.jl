###AGAL
#input
# u           vector to multiply
# xy          Q2 nodal coordinate vector
# xyp         Q1 nodal coordinate vector
# mv          Q2 element mapping matrix
# bound       indices of boundary points
# ae          local Q2 diffusion derivative matrix
#output
# w           A * u
function ho_agal(u, params)
    mesh = params["mesh"]; order = params["order"]
    dof = params["dof"]; ne = params["ne"]
    bdy = params["bdy"]
    refel = params["refel"]
    centers = params["centers"]
    w = zeros(length(u))

    # zero dirichlet bdy conditions
    uu = copy(u)
    u[bdy] = zeros(length(bdy))

    # loop over elements
    for e in 1:ne
        idx = Mesh.get_node_indices(mesh, e, order)
        pts = Mesh.element_nodes(mesh, e, refel)
        (detJac, Jac) = Mesh.geometric_factors(mesh, refel, pts)
        brinkman_pts = Mesh.brinkman_tensor(pts, centers)
        eMat = Mesh.element_stiffness_brinkman(mesh, e, refel, detJac, Jac, brinkman_pts)
        w[idx] += eMat * u[idx]
    end

    w[bdy] = uu[bdy]
    return vec(w[1:dof])
end
