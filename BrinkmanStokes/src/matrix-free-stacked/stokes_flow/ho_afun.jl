###HO_AFUN matrix-free A block operator for high order stokes
function ho_afun(u, params)
    u = share(u)
    mesh = params["mesh"]
    order = params["order"]
    dof = params["dof"]
    ne = params["ne"]
    bdy = params["bdy"]
    refel = params["refel"]
    centers = params["centers"]
    mv = share(params["mv"])
    w = zeros(length(u))

    # zero dirichlet bdy conditions
    uu = copy(u)
    u[bdy] = zeros(length(bdy))
    u[bdy + dof] = zeros(length(bdy))

    # loop over elements
    for e in 1:ne
        idx = Mesh.get_node_indices(mesh, e, order)
        # idx = vec(mv[e,:]')
        pts = Mesh.element_nodes(mesh, e, refel)
        (detJac, Jac) = Mesh.geometric_factors(mesh, refel, pts)
        brinkman_pts = Mesh.brinkman_tensor(pts, centers)
        eMat = Mesh.element_stiffness_brinkman(mesh, e, refel, detJac, Jac, brinkman_pts)
        w[idx] += eMat * u[idx]
        w[idx + dof] += eMat * u[idx + dof]
    end

    # @sync begin
    #   for p in procs()
    #     @async remotecall_wait(p, ho_afun_loop_chunk!, w, u, mesh, order, dof, refel, centers, mv)
    #   end
    # end

    w[bdy] = uu[bdy]
    w[bdy + dof] = uu[bdy + dof]
    return vec(w)
end

@everywhere function ho_afun_loop!(w, u, mesh, order, dof, refel, centers, mv, prange)
    pnel = length(prange)
    mve = mv[prange, :]
    for e in 1:pnel
        # idx = Mesh.get_node_indices(mesh, e, order)
        idx = vec(mve[e, :]')
        pts = Mesh.element_nodes(mesh, e, refel)
        detJac, Jac = Mesh.geometric_factors(mesh, refel, pts)
        brinkman_pts = Mesh.brinkman_tensor(pts, centers)
        eMat = Mesh.element_stiffness_brinkman(mesh, e, refel, detJac, Jac, brinkman_pts)
        w[idx] += eMat * u[idx]
        w[idx + dof] += eMat * u[idx + dof]
    end
    w
end

@everywhere ho_afun_loop_chunk!(w, u, mesh, order, dof, refel, centers, mv) = ho_afun_loop!(w, u, mesh, order, dof, refel, centers, mv, myrange(mv))
