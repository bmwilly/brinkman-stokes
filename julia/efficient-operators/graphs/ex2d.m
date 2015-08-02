% specify the coefficients
mu = @(x,y)(1 + 1e6*( cos(2*pi*x)^2 + cos(2*pi*y)^2 ) );

% create the mesh heirarchy, for a warped mesh
% in this case we create a h+p heirarchy
%     grid 4 --> 16x16, p=4  (finest)
%     grid 3 --> 16x16, p=2      
%     grid 2 --> 16x16, p=1
%     grid 1 -->  8x8 , p=1  (coarsest)

disp('--------------------');
disp('Creating Mesh heirarchy');
% g = create_hexmesh_grids(2, mu, @homg.xform.shell, [1 2 4], [8 16]);

dim = 2; xform = @homg.xform.identity; orders = [1 2 4]; nelems = [8 16];

num_hgrids = length(nelems);
num_pgrids = length(orders);

num_grids = num_hgrids + num_pgrids - 1;

m = homg.hexmesh(repmat(nelems(1), 1, dim), xform);
coarse = homg.grid(m, orders(1));

% create h-grids
for i=2:num_hgrids
  m = homg.hexmesh(repmat(nelems(i), 1, dim), xform);
  grid = homg.grid(m, orders(1), coarse);
  coarse = grid;
end

hfine = nelems(num_hgrids);

% create p-grids
for i=2:num_pgrids
  m = homg.hexmesh(repmat(hfine, 1, dim), xform);
  grid = homg.grid(m, orders(i), coarse);
  coarse = grid;
end

% Now assemble matrices ...
grid.assemble_poisson(mu);

% g.solve_pcg(# vcycles, smoother, pre-smoothing steps, post-smoothing steps, rhs, initial guess)

% disp('--------------------');
% disp('==== Solving using Multigrid as a solver ====');
% % now solve using multigrid as the solver and the choice of smoother
% disp('==== Jacobi ====');
% g.solve (150, 'jacobi', 3,3, g.L, g.get_u0);
% disp('==== Chebyshev ====');
% g.solve (150, 'chebyshev', 3,3, g.L, g.get_u0);
% disp('==== SSOR ====');
% g.solve (150, 'ssor', 2,1, g.L, g.get_u0);

% disp('--------------------');
% disp('==== Solving using Multigrid as a preconditioner ====');
% % or solve using CG preconditioned using multigrid
% disp('==== Jacobi ====');
% g.solve_pcg(150, 'jacobi', 3,3, g.L, g.get_u0);
% disp('==== Chebyshev ====');
% g.solve_pcg(150, 'chebyshev', 3,3, g.L, g.get_u0);
% disp('==== SSOR ====');
% g.solve_pcg(150, 'ssor', 2,1, g.L, g.get_u0);
% disp('--------------------');

