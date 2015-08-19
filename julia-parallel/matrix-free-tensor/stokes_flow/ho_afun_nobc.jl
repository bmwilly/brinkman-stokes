###HO_AFUN matrix-free A block operator for high order stokes
function ho_afun_nobc(u, params)
  
  mesh = params["mesh"]; order = params["order"];
  dof = params["dof"]; ne = params["ne"]; NP = params["NP"];
  bdy = params["bdy"]; eMat = params["eMat"];
  w = zeros(length(u))

  Ux = zeros(NP, ne); Uy = zeros(NP, ne)

  # loop over elements
  for e = 1:ne
    idx = Mesh.get_node_indices(mesh, e, order)
    Ux[:, e] = u[idx]
    Uy[:, e] = u[idx+dof]
  end
  Wx = eMat * Ux; Wy = eMat * Uy
  for e = 1:ne
    idx = Mesh.get_node_indices(mesh, e, order)
    w[idx] += Wx[:, e]
    w[idx+dof] += Wy[:, e]
  end
  vec(w)
end
