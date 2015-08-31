function [a,b,bbx,bby,f,g] = stokes_q2p1(xy,xyp,mv)
%STOKES_Q2P1 Q2-P1 matrix generator
%   [A,B,Bx,By,f,g] = stokes_q2p1(xy,xyp,mv);
%   input
%          xy         Q2 nodal coordinate vector 
%          xyp        Q1 nodal coordinate vector  
%          mv         Q2 element mapping matrix
%   output
%          A          Q2 vector diffusion matrix
%          B          Q2-Q1 divergence matrix 
%          Bx         Q2 x-derivative matrix    
%          By         Q2 y-derivative matrix    
%          f          velocity rhs vector
%          g          pressure rhs vector
%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function flowbc.

nngpt=9;     
x=xy(:,1); 
y=xy(:,2);
xp=xyp(:,1); 
yp=xyp(:,2);
nvtx=length(x); 
nu=2*nvtx; 
np=3*length(xp); 
nel=length(mv(:,1));  
mp=[[1:3:3*nel]',[2:3:3*nel]',[3:3:3*nel]'];
lx=max(x)-min(x); 
ly=max(y)-min(y);
hx=max(diff(x)); 
hy=max(diff(y));
fprintf('setting up Q2-P1 matrices...  ')
%
% initialise global matrices
a = sparse(nu,nu);
bbx = sparse(nvtx,nvtx);
bby = sparse(nvtx,nvtx);
bx = sparse(np,nvtx);
by = sparse(np,nvtx);
b = sparse(np,nu);
f = zeros(nu,1);
g = zeros(np,1);
%
%
% Gauss point integration rules
% 3x3 Gauss points
gpt=sqrt(0.6); 
s(1) = -gpt; t(1) = -gpt; wt(1)=25/81;
s(2) =  gpt; t(2) = -gpt; wt(2)=25/81;
s(3) =  gpt; t(3) =  gpt; wt(3)=25/81; 
s(4) = -gpt; t(4) =  gpt; wt(4)=25/81;
s(5) =  0.0; t(5) = -gpt; wt(5)=40/81;
s(6) =  gpt; t(6) =  0.0; wt(6)=40/81;
s(7) =  0.0; t(7) =  gpt; wt(7)=40/81; 
s(8) = -gpt; t(8) =  0.0; wt(8)=40/81;
s(9) =  0.0; t(9) =  0.0; wt(9)=64/81;
%
% inner loop over elements    
for ivtx = 1:4
    xl_v(:,ivtx) = x(mv(:,ivtx));
    yl_v(:,ivtx) = y(mv(:,ivtx)); 
end

ae = zeros(nel,9,9);
bbxe = zeros(nel,9,9);
bbye = zeros(nel,9,9);
bxe = zeros(nel,3,9);
bye = zeros(nel,3,9);
ge = zeros(nel,3);
% 
% loop over Gauss points
for igpt = 1:nngpt
    sigpt=s(igpt);
    tigpt=t(igpt);
    wght=wt(igpt);
    
    %  evaluate derivatives etc
    [jac,invjac,phi,dphidx,dphidy] = deriv(sigpt,tigpt,xl_v,yl_v);
    [psi,dpsidx,dpsidy] = qderiv(sigpt,tigpt,xl_v,yl_v);        
    [chi,dchidx,dchidy] = lderiv(sigpt,tigpt,xl_v,yl_v);
    for j = 1:9
        for i = 1:9
            ae(:,i,j)  = ae(:,i,j)  + wght*dpsidx(:,i).*dpsidx(:,j).*invjac(:);
            ae(:,i,j)  = ae(:,i,j)  + wght*dpsidy(:,i).*dpsidy(:,j).*invjac(:);
            bbxe(:,i,j) = bbxe(:,i,j) - wght*psi(:,i) .*dpsidx(:,j);             
            bbye(:,i,j) = bbye(:,i,j) - wght*psi(:,i) .*dpsidy(:,j);   
        end
        for i=1:3
            bxe(:,i,j) = bxe(:,i,j) - wght*chi(:,i) .* dpsidx(:,j);
            bye(:,i,j) = bye(:,i,j) - wght*chi(:,i) .* dpsidy(:,j);
        end
    end
    %
% end of Gauss point loop
end  
%
%%  element assembly into global matrices
% component velocity matrices ...    
for krow=1:9
    nrow=mv(:,krow);	 
    for kcol=1:9
        ncol=mv(:,kcol);	  
        a = a + sparse(nrow, ncol, ae(:,krow,kcol), nu, nu);
        a = a + sparse(nrow+nvtx, ncol+nvtx, ae(:,krow,kcol), nu, nu);
        bbx = bbx + sparse(nrow, ncol, bbxe(:,krow,kcol), nvtx, nvtx);
        bby = bby + sparse(nrow, ncol, bbye(:,krow,kcol), nvtx, nvtx);
    end
    for kcol=1:3
        ncol=mp(:,kcol);	  
        bx = bx + sparse(ncol, nrow, bxe(:,kcol,krow), np, nvtx);
        by = by + sparse(ncol, nrow, bye(:,kcol,krow), np, nvtx);
    end
end

%
% vector velocity matrices ...
b = [bx,by];
%
fprintf('done\n')
return

