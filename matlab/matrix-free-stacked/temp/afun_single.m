function w = afun_single(u, p)

nu = p.nu; np = p.np; ae = p.ae; mv = p.mv; nel = p.nel; nvtx = p.nvtx;

aes = squeeze(ae(1,:,:));

w = zeros(nu+np,1);
[n,m] = size(mv);

% for e = 1:nel
%     ind = mv(e,:);
%     we = aes * u(ind);
%     w(ind) = w(ind) + we;
%     w(ind + nvtx) = w(ind + nvtx) + we;
% end

Ux = zeros(m,n);
Uy = zeros(m,n);
for e = 1:nel
    ind = mv(e,:);
    Ux(:,e) = u(ind);
    Uy(:,e) = u(ind+nvtx);
end

Wx = aes * Ux;
Wy = aes * Uy;

for e = 1:nel
    ind = mv(e,:);
    w(ind) = w(ind) + Wx(:,e);
    w(ind+nvtx) = w(ind+nvtx) + Wy(:,e);
end

end
