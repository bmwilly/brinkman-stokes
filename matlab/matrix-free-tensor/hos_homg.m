function etoc = hos_homg(order, msize, dim)
  nelems = [2^msize];
  m = homg.hexmesh(repmat(nelems, 1, dim), homg.xform.identity);
  refel = homg.refel(m.dim, order);
  dof = prod(m.nelems*order + 1);
  ne = prod(m.nelems);
  NP = (order+1)^m.dim;
  NPNP = NP*NP;
  bdy = m.get_boundary_node_indices(order);
  pts = m.element_nodes(1, refel);
  [detJac, Jac] = m.geometric_factors(refel, pts);
  eMat = m.element_stiffness(1, refel, detJac, Jac);

  params = {
    'mesh', m,
    'order', order,
    'dof', dof,
    'ne', ne,
    'NP', NP,
    'bdy', bdy,
    'eMat', eMat
  }

  tic;
  for cnt = 1:100; u = rand(2*dof, 1); w = ho_afun(u, params); end;
  etoc = toc;
end
