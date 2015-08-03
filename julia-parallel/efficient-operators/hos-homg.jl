reload("helpers/input.jl")
reload("../julia-homg/Basis.jl")
reload("../julia-homg/Hexmesh.jl")
reload("../julia-homg/Xform.jl")
reload("../julia-homg/Grids.jl")
reload("../julia-homg/Tensor.jl")
reload("../julia-homg/Refel.jl")

dim = 2
order = int(input("Polynomial order: "))
msize = int(input("Mesh size: "))
# msize = 2
nelems = [2^msize]

# nelems = [msize]; dim = 2; order = 2;

m = Mesh.Hexmesh(tuple(repmat(nelems, 1, dim)...), Xform.identity)
# Mesh.plot(m)
dof = prod([m.nelems...]*order + 1)
tic()
K,M,iK = Mesh.assemble_poisson(m, order)
k1,k2 = size(K)
A = [K spzeros(k1,k2); spzeros(k1,k2) K]
println(toc())
@time for cnt = 1:100; u = vec(rand(2dof, 1)); w = A*u; end

# u = linspace(1, 2dof, 2dof)
# w = A*u


# self = m;
#
# Mesh.set_order(self,order);
# # assemble the mass matrix
# refel = Refel( self.dim, order );
# dof = prod([self.nelems...]*order + 1);
# ne = prod([self.nelems...]);
# # storage for indices and values
# NP = (order+1)^self.dim;
# NPNP = NP * NP;
#
# I = zeros(ne * NPNP, 1);
# J = zeros(ne * NPNP, 1);
# mass_val = zeros(ne * NPNP, 1);
# stiff_val = zeros(ne * NPNP, 1);
# inv_stiff_val = zeros(ne * NPNP, 1);
# ind_inner1D = repmat((2:order), 1, order-1);
# if self.dim == 2
#
#   ind_inner = ind_inner1D + (order+1) * (ind_inner1D'-1);
# else
#   ind_inner = ind_inner1D + (order+1) * (ind_inner1D'-1);
#   ind_inner = repmat(ind_inner, [1,1,order-1]);
#   for i = 1:order-1
#     ind_inner[:,:,i] = ind_inner[:,:,i] + i * (order+1)^2;
#   end
# end
#
# # loop over elements
# for e=1:ne
#   idx =  Mesh.get_node_indices(self, e, order);
#   ind1 = repmat(idx,NP,1);
#   ind2 = reshape(repmat(idx',NP,1),NPNP,1);
#   st = (e-1)*NPNP+1;
#   en = e*NPNP;
#
#   I[st:en] = ind1;
#   J[st:en] = ind2;
#   pts =  Mesh.element_nodes(self, e, refel);
#   (detJac, Jac) = Mesh.geometric_factors(self, refel, pts);
#
#   eMat = Mesh.element_mass(self, e, refel, detJac);
#   mass_val[st:en] = eMat[:];
#
#   eMat = Mesh.element_stiffness(self, e, refel, detJac, Jac);
#   stiff_val[st:en] = eMat[:];
#
#   eMat_inner_inv = inv(eMat[ind_inner[:],ind_inner[:]]);
#   eMat_inv = diagm(diag(eMat,0));
#
#   eMat_inv[ind_inner[:],ind_inner[:]] =  eMat_inner_inv;
#   inv_stiff_val[st:en] = eMat_inv[:];
# end
#
# Iv=int64(I[:]);
# Jv=int64(J[:]);
# mv=mass_val[:];
#
# M = sparse(Iv,Jv,mv,dof,dof);
# # zero dirichlet bdy conditions
# bdy = Mesh.get_boundary_node_indices(self, order);
# ii = Mesh.ismember(I,bdy);
# jj = Mesh.ismember(J,bdy);
#
# stiff_val = stiff_val.*(int(!bool(ii))).*(int(!bool(jj)));
# inv_stiff_val = inv_stiff_val.*(int(!bool(ii))).*(int(!bool(jj)));
#
# I = [I; bdy];
# J = [J; bdy];
# stiff_val = [stiff_val; ones(length(bdy), 1)];
# inv_stiff_val = [inv_stiff_val; ones(length(bdy), 1)];
# Iv=int64(I[:]);
# Jv=int64(J[:]);
# sv=stiff_val[:];
# isv=inv_stiff_val[:];
#
# K = sparse(Iv,Jv,sv,dof,dof);
# iK = sparse(Iv,Jv,isv,dof,dof);
# ebdy = Mesh.get_element_boundary_node_indices(self, order);
# iKebdry = diag(full(iK[ebdy,ebdy]),0)
# if countnz(iKebdry) > 0
#     iK[ebdy,ebdy] = diagm(1./iKebdry)
# end
# return K, M, iK
#
# g = Grids.Grid(m, order)
# mu = (x,y) -> 1
# Grids.assemble_poisson(g, mu)
