###HO_AFUN matrix-free A block operator for high order stokes
function ho_afun(u, params)
  u = share(u)
  mesh = params["mesh"]; order = params["order"];
  dof = params["dof"]; ne = params["ne"]; NP = params["NP"];
  bdy = params["bdy"];
  refel = params["refel"]
  centers = params["centers"]
  mv = share(params["mv"])
  w = zeros(length(u))
  Ux = zeros(NP, ne); Uy = zeros(NP, ne)
  eMats = zeros(NP*ne, NP)

  # zero dirichlet bdy conditions
  uu = copy(u)
  u[bdy] = zeros(length(bdy))
  u[bdy+dof] = zeros(length(bdy))

  # loop over elements
  for e = 1:ne
    idx = Mesh.get_node_indices(mesh, e, order)
    Ux[:,e] = u[idx]
    Uy[:,e] = u[idx+dof]
    pts = Mesh.element_nodes(mesh, e, refel)
    (detJac, Jac) = Mesh.geometric_factors(mesh, refel, pts)
    brinkman_pts = Mesh.brinkman_tensor(pts, centers);
    eMat = Mesh.element_stiffness_brinkman(mesh, e, refel, detJac, Jac, brinkman_pts)
    eMats[(e-1)*NP+1:e*NP,:] = eMat
  end
  Wx = eMats * Ux; Wy = eMats * Uy;
  for e = 1:ne
    idx = Mesh.get_node_indices(mesh, e, order)
    w[idx] += Wx[(e-1)*NP+1:e*NP,e]
    w[idx+dof] += Wy[(e-1)*NP+1:e*NP,e]
  end

  w[bdy] = uu[bdy]
  w[bdy+dof] = uu[bdy+dof]
  vec(w)
end
