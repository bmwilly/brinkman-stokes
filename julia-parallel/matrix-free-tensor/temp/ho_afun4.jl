channel_grid = channel_domain(msize) # Q2 grid for channel domain
stokes_grid = q2p1grid(channel_grid)
stokes_mats = stokes_q2p1(stokes_grid) # stokes element matrices

bounds = {
  "bound" => channel_grid["bound"],
  "bndxy" => channel_grid["bndxy"],
  "bnde" => channel_grid["bnde"],
  "obs" => channel_grid["obs"]
}

# brinkman obstacles
centers = [0.33 0.5]

order = 2; dim = 2;
nelems = [2^(msize)]
m = Mesh.Hexmesh(tuple(repeat(nelems, 1, dim)...), Xform.identity)
dof = prod([m.nelems...] * order + 1)
Mesh.set_order(m,order);
refel = Refel(m.dim, order);
dof = prod([m.nelems...] * order + 1);
ne = prod([m.nelems...]);
# storage for indices and values
NP = (order + 1)^m.dim;
NPNP = NP * NP;
I = zeros(ne * NPNP, 1);
J = zeros(ne * NPNP, 1);
stiff_val = zeros(ne * NPNP, 1);
bdy = Mesh.get_boundary_node_indices(m, order);
# pts =  Mesh.element_nodes(m, 1, refel);
# (detJac, Jac) = Mesh.geometric_factors(m, refel, pts);
# brinkman_pts = Mesh.brinkman_tensor(pts, centers)

w = zeros(length(u))

# zero dirichlet bdy conditions
u = linspace(1, 2dof, 2dof);
uu = copy(u)
u[bdy] = zeros(length(bdy))
u[bdy + dof] = zeros(length(bdy))
Ux = zeros(NP, ne); Uy = zeros(NP, ne)

# loop over elements
for e = 1:ne
    idx = Mesh.get_node_indices(mesh, e, order)
    ind1 = repeat(idx, NP, 1)
    ind2 = reshape(repeat(idx', NP, 1), NPNP, 1);
    st = (e - 1) * NPNP + 1;
    en = e * NPNP;
    I[st:en] = ind1;
    J[st:en] = ind2;
    pts = Mesh.element_nodes(mesh, e, refel)
    (detJac, Jac) = Mesh.geometric_factors(mesh, refel, pts)
    brinkman_pts = Mesh.brinkman_tensor(pts, centers);
    eMat = Mesh.element_stiffness_brinkman(mesh, e, refel, detJac, Jac, brinkman_pts)
    stiff_val[st:en] = eMat[:]

  # w[st:en] += eMat * u[st:en]
    w[idx] += eMat * u[idx]

#   Ux[:, e] = u[idx]
#   Uy[:, e] = u[idx+dof]
# end
# Wx = eMat * Ux; Wy = eMat * Uy
# for e = 1:ne
#   idx = Mesh.get_node_indices(mesh, e, order)
#   w[idx] += Wx[:, e]
#   w[idx+dof] += Wy[:, e]
end

# ii = ismember(I,bdy);
# jj = ismember(J,bdy);
# stiff_val = stiff_val.*(int(!bool(ii))).*(int(!bool(jj)));
# I = [I; bdy];
# J = [J; bdy];
# stiff_val = [stiff_val; ones(length(bdy), 1)];
# Iv=int64(I[:]);
# Jv=int64(J[:]);
# sv=stiff_val[:];
# K = sparse(Iv,Jv,sv,dof,dof);

# for e = 1:ne
#   idx = Mesh.get_node_indices(mesh, e, order)
#   ind1 = repeat(idx,NP,1)
#   ind2 = reshape(repeat(idx',NP,1),NPNP,1);
#   st = (e-1)*NPNP+1;
#   en = e*NPNP;
#   I[st:en] = ind1;
#   J[st:en] = ind2;
#   eMat = Mesh.element_mass_brinkman(mesh, e, refel, detJac, brinkman_pts)
#   Wx = eMat * Ux
#   w[idx] += Wx[:,e]
# end

w[bdy] = uu[bdy]
w[bdy + dof] = uu[bdy + dof]
vec(w)
