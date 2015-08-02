function f = streambc(f, xy, xyp, mv, bound, ae)
%STREAMBC imposes Dirichlet BC on the streamfunction

nvtx = length(f); 
fx=f(1:nvtx);

%% set boundary condition
xbd=xy(bound,1); 
bcx=0*xbd;
fx(bound)=bcx; 
f=fx;

end
