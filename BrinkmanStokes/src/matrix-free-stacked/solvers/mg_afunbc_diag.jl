###MG_AFUNBC_DIAG modified diagonal operator with zero boundary condition imposed
function mg_afunbc_diag(u, mparams)
  # mesh = params["mesh"];
  # order = params["order"];
  # dof = params["dof"];
  # ne = params["ne"];
  # bdy = params["bdy"];
  # refel = params["refel"];
  # centers = params["centers"];
  # ev = share(params["ev"])
  # w = zeros(length(u))

  ev = mparams["ev"]; bound = mparams["bound"];
  ae = mparams["ae"]
  nel = length(ev[:, 1])
  aes = squeeze(ae[1, :, :], 1)
  w = zeros(length(u))

  # zero dirichlet boundary conditions
  uu = copy(u)
  u[bound] = zeros(length(bound))

  # for e = 1:ne
  #   idx = Mesh.get_node_indices(mesh, e, order)
  #   pts = Mesh.element_nodes(mesh, e, refel)
  #   brinkman_pts = Mesh.brinkman_tensor(pts, centers);
  #   eMat = Mesh.element_stiffness_brinkman(mesh, e, refel, detJac, Jac, brinkman_pts)
  #   w[idx] += diagm(diag(MM)) * u[idx]
  # end

  for e = 1:nel
    ind = vec(ev[e, :]')
    w[ind] += diagm(diag(aes)) * u[ind]
  end

  w[bound] = uu[bound]
  vec(w)
end
