function w = afunbc(u, xy, xyp, mv, bound, ae)

x = xy(:, 1);
xp = xyp(:, 1);
nvtx = length(x);
nu = 2 * nvtx;
np = 3 * length(xp);
nel = length(mv(:, 1));
aes = squeeze(ae(1,:,:));

% zero dirichlet boundary condition
uu = u;
u(bound) = zeros(length(bound), 1);
u(bound+nvtx) = zeros(length(bound), 1); 

w = zeros(nu+np,1);
for e = 1:nel
    ind = mv(e,:);
    w(ind) = w(ind) + aes * u(ind); 
    w(ind+nvtx) = w(ind+nvtx) + aes * u(ind+nvtx);
end

w(bound) = uu(bound); 
w(bound+nvtx) = uu(bound+nvtx);

end
