function w = bbfunplot(u, xy, xyp, mv, bbxe, bbye)

x = xy(:, 1);
xp = xyp(:, 1);
nvtx = length(x);
nu = 2 * nvtx;
np = 3 * length(xp);
nel = length(mv(:, 1));

ux = u(1:nvtx); 
uy = u(nvtx + 1:nu);

bbxes = squeeze(bbxe(1, :, :));
bbyes = squeeze(bbye(1, :, :));

w = zeros(nvtx, 1); 
for e = 1:nel
    ind = mv(e, :);
    wex = bbxes * uy(ind);
    wey = bbyes * ux(ind);
    w(ind) = w(ind) + wey;
    w(ind) = w(ind) - wex;
end

end
