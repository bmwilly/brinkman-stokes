function w = btfun(u, xy, xyp, mv, bxe, bye)

x = xy(:, 1);
xp = xyp(:, 1);
nvtx = length(x);
nu = 2 * nvtx;
np = 3 * length(xp);
nel = length(mv(:, 1));
mp =[[1:3:3*nel]',[2:3:3*nel]',[3:3:3*nel]'];

up = u(nu+1:nu+np);

bxes = squeeze(bxe(1,:,:));
byes = squeeze(bye(1,:,:));

w = zeros(nu+np,1);
for e = 1:nel
    ind9 = mv(e,:); 
    ind3 = mp(e,:);
    
    wepx = up(ind3)' * bxes; 
    wepy = up(ind3)' * byes; 
    
    w(ind9) = w(ind9) + wepx';
    w(ind9+nvtx) = w(ind9+nvtx) + wepy';
end

end