###HO_AFUN matrix-free A block operator for high order stokes
function ho_afun_square(u, params)
  mesh = params["mesh"];
  order = params["order"];
  dof = params["dof"];
  ne = params["ne"];
  bdy = params["bdy"];
  refel = params["refel"];
  mv = share(params["mv"])
  w = zeros(length(u))

  # zero dirichlet bdy conditions
  uu = copy(u)
  u[bdy] = zeros(length(bdy))
  u[bdy+dof] = zeros(length(bdy))

  # loop over elements
  for e = 1:ne
    idx = Mesh.get_node_indices(mesh, e, order)
    pts = Mesh.element_nodes(mesh, e, refel)
    (detJac, Jac) = Mesh.geometric_factors(mesh, refel, pts)
    eMat = Mesh.element_stiffness(mesh, e, refel, detJac, Jac)
    w[idx] += eMat * u[idx]
    w[idx+dof] += eMat * u[idx+dof]
  end

  w[bdy] = uu[bdy]
  w[bdy+dof] = uu[bdy+dof]
  vec(w)
end
