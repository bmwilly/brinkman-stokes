function t = mv_fun()
%     solve_square_stokes
    square_stokes
    stokes_qp
    tic;
    for cnt = 1:100
      u = rand(nu + np, 1);
%       w = kfunbc(u, xy, xyp, mv, bound, ae, bxe, bye);
      w = afunbc(u, xy, xyp, mv, bound, ae);
    end
    t = toc;
end