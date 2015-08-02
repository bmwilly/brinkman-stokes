%SOLVE_STOKES solve Stokes problem in square domain
%
fprintf('imposing (enclosed flow) boundary conditions ...\n')
tic;
%% compute solution
%-------------------------------------------------------
% matrix-free method
fst = fbc(xy, xyp, mv, bound, ae, bxe, bye);
gst = gbc(xy, xyp, mv, bound, bxe, bye);
xst = minres(@(u)kfunbc(u, xy, xyp, mv, bound, ae, bxe, bye), [fst; gst], 1e-6, max(length([fst; gst]), 1000));
%-------------------------------------------------------
etoc=toc;
fprintf('Stokes system solved in %8.3e seconds\n',etoc)
%
