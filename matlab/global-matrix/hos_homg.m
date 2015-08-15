function etoc = hos_homg(order, msize, dim)
  nelems = [2^msize];
  m = homg.hexmesh(repmat(nelems, 1, dim), homg.xform.identity);
  dof = prod(m.nelems*order + 1);
  K = homg.assemble_poisson(m, order);
  [k1, k2] = size(K);
  A = [K spzeros(k1,k2); spzeros(k1,k2) K];
  tic;
  for cnt = 1:100; u = rand(2*dof); w = A*u; end;
  etoc = toc;
end
