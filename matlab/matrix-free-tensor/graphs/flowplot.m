function flowplot(sol, xy, xyp, mv, bound, ae, bxe, bye, bbxe, bbye)
%FLOWPLOT 
%   flowplot(sol,By,Bx,A,xy,xyp,x,y,bound);
%   input
%          sol        flow solution vector
%          By         velocity  y-derivative matrix    
%          Bx         velocity x-derivative matrix    
%          A          vector diffusion matrix
%          xy         velocity nodal coordinate vector  
%          xyp        pressure nodal coordinate vector  
%          x          vector of x-axis interpolation points
%          y          vector of y-axis interpolation points
%          bound      boundary vertex vector
%
% calls function streambc.m to set boundary values

x = xy(:, 1);
y = xy(:, 2);
xp = xyp(:, 1);
yp = xyp(:, 2);
nvtx=length(xy); 
nu=2*nvtx; 
np=length(xyp);
%
% compute auxilliary quantites
u=sol(1:nu);
p=sol(nu+1:end);
f = bbfunplot(u, xy, xyp, mv, bbxe, bbye);
fsv = streambc(f, xy, xyp, mv, bound, ae);
phi = minres(@(u)afunstreambc(u, xy, xyp, mv, bound, ae), fsv, 1e-9, 500); 
%
%% plot pressure
p=p(1:3:end); 
xx=x(1:end); 
yy=y(1:end);

% interpolate to a cartesian product mesh
[X,Y]=meshgrid(xx,yy);
xysol = griddata(xyp(:,1),xyp(:,2),p,X,Y);
figure
subplot(122), mesh(X,Y,xysol),axis('square')
title('pressure field','FontSize',12)
%
%% plot velocity
ax = [min(x)-.1 max(x)+.1 min(y)-.1 max(y)+.1];
[X,Y]=meshgrid(x,y);
xysol = griddata(xy(:,1),xy(:,2),phi,X,Y);
maxphi=max(max(xysol)); 
minphi=min(min(xysol));
subplot(121),contour(X,Y,xysol,24),axis('square'), axis(ax)	
title('Streamlines: uniform','FontSize',12); 
hold on
plot([-1,1],[-1,-1],'-k');
plot([1,1],[-1,1],'-k');
plot([1,-1],[1,1],'-k');
plot([-1,-1],[1,-1],'-k');
hold off
axis('off')
return
