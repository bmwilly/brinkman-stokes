%SOLVE_STOKES solve Stokes problem in square domain
%
fprintf('imposing (enclosed flow) boundary conditions ...\n')
%% boundary conditions
[Ast,Bst,fst,gst] = flowbc(A,B,f,g,xy,bound);
%
np=length(gst);
nu=length(fst);
tic;
%% compute solution
%-------------------------------------------------------
K = [Ast, Bst'; Bst, sparse(np,np)];
[kn, km] = size(K);
xst = minres(K, [fst; gst], 1e-6, max(kn, 1000));
%-------------------------------------------------------
etoc=toc;
fprintf('Stokes system solved in %8.3e seconds\n',etoc)
%
%%% plot solution
% flowplot(xst,By,Bx,A,xy,xyp,x,y,bound);
fprintf('\n')
