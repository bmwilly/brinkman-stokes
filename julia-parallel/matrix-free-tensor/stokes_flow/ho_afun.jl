###HO_AFUN matrix-free A block operator for high order stokes
function ho_afun(u, params)

  mesh = params["mesh"]; order = params["order"];
  dof = params["dof"]; ne = params["ne"]; NP = params["NP"];
  bdy = params["bdy"]; eMat = params["eMat"];
  # w = zeros(length(u))
  # w = SharedArray(Float64, length(u))
  eMat = share(eMat)
  w = share(zeros(length(u)))

  # zero dirichlet bdy conditions
  uu = copy(u)
  u[bdy] = zeros(length(bdy))
  u[bdy+dof] = zeros(length(bdy))

  # Ux = zeros(NP, ne); Uy = zeros(NP, ne)
  # Ux = SharedArray(Float64, (NP, ne))
  # Uy = SharedArray(Float64, (NP, ne))
  Ux = share(zeros(NP, ne)); Uy = share(zeros(NP, ne));

  # loop over elements
  # @sync @parallel for e = 1:ne
  for e = 1:ne
    idx = Mesh.get_node_indices(mesh, e, order)
    Ux[:, e] = u[idx]
    Uy[:, e] = u[idx+dof]
  end
  Wx = eMat * Ux; Wy = eMat * Uy
  # @sync @parallel for e = 1:ne
  for e = 1:ne
    idx = Mesh.get_node_indices(mesh, e, order)
    w[idx] += Wx[:, e]
    w[idx+dof] += Wy[:, e]
  end

  w[bdy] = uu[bdy]
  w[bdy+dof] = uu[bdy+dof]
  vec(w)
end
