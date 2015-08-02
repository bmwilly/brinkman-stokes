function w = gbc(xy, xyp, mv, bound, bxe, bye)  

% get variables
x = xy(:, 1);
xp = xyp(:, 1);
nvtx = length(x);
nu = 2 * nvtx;
np = 3 * length(xp);

% create vector without imposed boundary conditions
% w = bfun(u, p);
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
w = w - bfun(bc, xy, xyp, mv, bxe, bye);

w = w(nu+1:nu+np);

end