function [x,y,xy,xyp] = q2p1grid(x,y,xy,mv,bound);
%Q2P1GRID Q2-P1 element grid generator
%   [x,y,xy,xyp] = q2p1grid(x,y,xy,mv,bound);
%   input
%          x          x coordinate vector
%          y          y coordinate vector 
%          xy         nodal coordinate vector  
%          mv         Q2 macroelement mapping matrix
%          bound      boundary vertex vector
%   output       
%          xyp        centroid coordinate vector
xx=xy(:,1); 
yy=xy(:,2); 
nvtx=length(xx);
nel=length(mv(:,1));

%% recompute mid-side points in the case of stretched grids 
% y-direction
yv=yy; 
ny=length(y);

for k=2:2:ny;
    yold=y(k); 
    ynew=0.5*(y(k+1)+y(k-1));
    l=find(yy==yold); 
    yv(l)=ynew; 
    y(k)=ynew;
end

% x-direction
xv=xx; 
nx=length(x);
for k=2:2:nx;
    xold=x(k); 
    xnew=0.5*(x(k+1)+x(k-1));
    l=find(xx==xold); 
    xv(l)=xnew; 
    x(k)=xnew;
end
xy=[xv,yv];

% centroid coordinates
for ielem=1:nel
    xc(ielem)=mean(xx(mv(ielem,1:4))); 
    yc(ielem)=mean(yy(mv(ielem,1:4)));
end
xyp=[xc',yc'];
%

% plotting of the grid 
% adj=sparse(nvtx,nvtx);
% for i=1:nel
%     adj(mv(i,1),mv(i,2))=1;
%     adj(mv(i,2),mv(i,3))=1;
%     adj(mv(i,3),mv(i,4))=1;
%     adj(mv(i,4),mv(i,1))=1;
% end
% figure(30)
% gplot(adj,xy,'b');
% axis('square')
% hold on
% plot(xy(:,1),xy(:,2),'ro')
% xybd=xy(bound,:);
% plot(xybd(:,1),xybd(:,2),'ko')
% plot(xyp(:,1),xyp(:,2),'k*',xyp(:,1),xyp(:,2),'ro')
% hold off
% title('Q2-P1 finite element subdivision')

return
