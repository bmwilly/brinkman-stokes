function w = afun(u, xy, xyp, mv, ae)

x = xy(:, 1);
xp = xyp(:, 1);
nvtx = length(x);
nu = 2 * nvtx;
np = 3 * length(xp);
nel = length(mv(:, 1));

aes = squeeze(ae(1,:,:));

w = zeros(nu+np,1);

for e = 1:nel
    ind = mv(e,:);
    we = aes * u(ind);
    w(ind) = w(ind) + we;
    w(ind + nvtx) = w(ind + nvtx) + we;
end

end
