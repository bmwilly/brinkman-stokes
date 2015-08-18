function etoc = hos_homg(order, msize, dim)
  nelems = [2^msize];
  m = homg.hexmesh(repmat(nelems, 1, dim), @homg.xform.identity);
  m.set_order(order);
  dof = prod(m.nelems*order + 1);
  K = m.assemble_poisson(order);
  [k1, k2] = size(K);
  A = [K sparse(k1,k2); sparse(k1,k2) K];
  tic;
  for cnt = 1:100; u = rand(2*dof,1); w = A*u; end;
  etoc = toc;
end
