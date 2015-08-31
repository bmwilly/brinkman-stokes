function w = fbc(xy, xyp, mv, bound, ae, bxe, bye)

% get variables
x = xy(:, 1);
xp = xyp(:, 1);
nvtx = length(x);
nu = 2 * nvtx;
np = 3 * length(xp);
nel = length(mv(:, 1));

% w = afun(u, p) + btfun(u, p);

% initialize vector
w = zeros(nu+np, 1);

% set boundary conditions 
xbd = xy(bound,1); 
ybd = xy(bound,2); 

% a regularized cavity
bcx=0*xbd;
bcy=0*xbd;
k=find(ybd==1 & xbd>-1 & xbd<1);
bcx(k)=(1-xbd(k).*xbd(k)).*(1+xbd(k).*xbd(k));

% impose boundary conditions
bccx = zeros(nvtx,1);
bccx(bound) = bcx;
bccy = zeros(nvtx,1);
bccy(bound) = bcy;

bc = [bccx; bccy; zeros(np,1)];
w = w - afun(bc, xy, xyp, mv, ae);
w = w - btfun(bc, xy, xyp, mv, bxe, bye);

wx = w(1:nvtx); 
wy = w(nvtx+1:nu); 
wp = w(nu+1:nu+np);

wx(bound) = bcx;
wy(bound) = bcy;

w = [wx;wy];

end