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
[n,m] = size(mv); 
Ux = zeros(m,n); Uy = zeros(m,n);

for e = 1:nel
    ind = mv(e,:);
    Ux(:,e) = u(ind); 
    Uy(:,e) = u(ind+nvtx);
end
Wx = aes * Ux; Wy = aes * Uy; 
for e = 1:nel
    ind = mv(e,:);
    w(ind) = w(ind) + Wx(:,e);
    w(ind+nvtx) = w(ind+nvtx) + Wy(:,e);
end

w(bound) = uu(bound); 
w(bound+nvtx) = uu(bound+nvtx);

end
