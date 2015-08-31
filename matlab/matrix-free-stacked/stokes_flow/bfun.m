function w = bfun(u, xy, xyp, mv, bxe, bye)

x = xy(:, 1);
xp = xyp(:, 1);
nvtx = length(x);
nu = 2 * nvtx;
np = 3 * length(xp);
nel = length(mv(:, 1));
mp = [[1:3:3*nel]',[2:3:3*nel]',[3:3:3*nel]'];

ux = u(1:nu/2);
uy = u(nu/2+1:nu);

bxes = squeeze(bxe(1,:,:));
byes = squeeze(bye(1,:,:));

w = zeros(nu+np,1);

for e = 1:nel
    ind9 = mv(e,:);
    wex = bxes * ux(ind9); 
    wey = byes * uy(ind9); 
    
    ind3 = mp(e,:) + nu; 
    w(ind3) = w(ind3) + wex; 
    w(ind3) = w(ind3) + wey;
end

end