addpath(genpath('/Users/bwilliams/Documents/brinkman-stokes/matlab/matrix-free-tensor'));
addpath(genpath('/Users/bwilliams/Documents/masters-thesis/MATLAB/homg'));

% order = 2; msize = 4; dim = 2; 
% t = hos_homg(order, msize, dim);

% orders = [2,3,4,5,6,7,8,9,10];
% orders = [11,12,13,14,15,16,17,18,19,20];
% orders = [21,22,23,24];
% msizes = [4];
msizes = [2,3,4,5,6,7,8,9,10];

% otimes = [];
% % dim = input('Dimension: ');
% dim = 2;
% if dim == 2
%   msize = 5;
% elseif dim == 3
%   msize = 3;
% end
% for order = orders
%   t = hos_homg(order, msize, dim);
% %   nelems = [2^msize];
% %   m = homg.hexmesh(repmat(nelems, 1, dim), @homg.xform.identity);
% %   refel = homg.refel(m.dim, order);
% %   dof = prod(m.nelems*order + 1);
% %   ne = prod(m.nelems);
% %   NP = (order+1)^m.dim;
% %   NPNP = NP*NP;
% %   bdy = m.get_boundary_node_indices(order);
% %   pts = m.element_nodes(1, refel);
% %   [detJac, Jac] = m.geometric_factors(refel, pts);
% %   eMat = m.element_stiffness(1, refel, detJac, Jac);
% % 
% %   params = struct('mesh', m, 'order', order,'dof', dof,'ne', ne, 'NP', NP,'bdy', bdy,'eMat', eMat)
% % 
% %   tic;
% %   for cnt = 1:100; u = rand(2*dof, 1); w = ho_afun(u, params); end;
% %   etoc = toc;
%   otimes = [otimes t];
% end

mtimes = []; 
dim = 2; 
order = 2;
for msize = msizes 
    t = hos_homg(order, msize, dim); 
    mtimes = [mtimes t];
end
disp mtimes
