function w = ho_afun(u, params)
  mesh = params.mesh; order = params.order;
  dof = params.dof; ne = params.ne; NP = params.NP;
  bdy = params.bdy; eMat = params.eMat;
  w = zeros(length(u), 1);

  uu = u;
  u(bdy) = zeros(length(bdy), 1);
  u(bdy+dof) = zeros(length(bdy), 1);

  Ux = zeros(NP, ne); Uy = zeros(NP, ne);

  for e = 1:ne
    idx = mesh.get_node_indices(e, order);
    Ux(:, e) = u(idx);
    Uy(:, e) = u(idx+dof);
  end
  Wx = eMat * Ux; Wy = eMat * Uy;
  for e = 1:ne
    idx = mesh.get_node_indices(e, order);
    w(idx) = w(idx) + Wx(:, e);
    w(idx+dof) = w(idx+dof) + Wy(:, e);
  end

  w(bdy) = uu(bdy);
  w(bdy+dof) = uu(bdy+dof);
end
