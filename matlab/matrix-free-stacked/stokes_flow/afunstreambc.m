function w = afunstreambc(u, xy, xyp, mv, bound, ae)

x = xy(:, 1);
nvtx = length(x);
nel = length(mv(:, 1));

aes = squeeze(ae(1, :, :));

w = zeros(nvtx, 1);
for e = 1:nel
    ind = mv(e, :);
    indbd = ismember(ind, bound);
    indint = ~ismember(ind, bound);
    indb = ind(indbd);

    ux = u(ind);
    ux(indbd) = zeros(1, length(find(indbd)));
    wex = aes * ux;

    % impose boundary conditions
    wex(indbd) = u(indb);
    
    % add interior points
    w(ind(indint)) = w(ind(indint)) + wex(indint);

    % add boundary points only when entry == 0
    wb = w(indb);
    ind0 = find(wb == 0);
    wexb = wex(indbd);
    w(indb(ind0)) = w(indb(ind0)) + wexb(ind0);
end

end