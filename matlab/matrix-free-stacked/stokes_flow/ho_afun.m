function w = ho_afun(u, params)
  mesh = params.mesh; order = params.order;
  dof = params.dof; ne = params.ne; NP = params.NP;
  bdy = params.bdy; eMat = params.eMat;
  w = zeros(length(u), 1);

  uu = u;
  u(bdy) = zeros(length(bdy), 1);
  u(bdy+dof) = zeros(length(bdy), 1);

  for e = 1:ne
    idx = m.get_node_indices(e, order)
    w(idx) = w(idx) + eMat * u(idx);
    w(idx+dof) = w(idx+dof) + eMat * u(idx+dof);
  end

  w(bdy) = uu(bdy);
  w(bdy+dof) = uu(bdy+dof);
end
