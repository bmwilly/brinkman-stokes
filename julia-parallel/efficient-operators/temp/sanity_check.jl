include("../../efficient-operators/stokes_flow/square_stokes.jl")
include("../../efficient-operators/stokes_flow/obstacle_stokes.jl")
include("../../efficient-operators/stokes_flow/brinkman_stokes.jl")
include("../../efficient-operators/stokes_flow/flowbc.jl")
include("../../efficient-operators/solvers/mg_diff.jl")
include("../../efficient-operators/solvers/m_st_mg.jl")
include("../hos_homg.jl")

dim = 2; msize = 5; order = 2;
nelems = [2^msize]

m = Mesh.Hexmesh(tuple(repeat(nelems, 1, dim)...), Xform.identity)
dof = prod([m.nelems...] * order + 1)

self = m;
Mesh.set_order(self,order);
# assemble the mass matrix
refel = Refel(self.dim, order);
dof = prod([self.nelems...] * order + 1);
ne = prod([self.nelems...]);
# storage for indices and values
NP = (order + 1)^self.dim;
NPNP = NP * NP;

I = zeros(ne * NPNP, 1);
J = zeros(ne * NPNP, 1);
mass_val = zeros(ne * NPNP, 1);
stiff_val = zeros(ne * NPNP, 1);
inv_stiff_val = zeros(ne * NPNP, 1);
ind_inner1D = repeat((2:order), 1, order - 1);
if self.dim == 2

    ind_inner = ind_inner1D + (order + 1) * (ind_inner1D' - 1);
else
    ind_inner = ind_inner1D + (order + 1) * (ind_inner1D' - 1);
  # ind_inner = repeat(ind_inner, [1,1,order-1]);
    ind_inner = repeat(ind_inner, outer=[1,1,order - 1]);
    for i = 1:order - 1
    # ind_inner[:,:,i] = ind_inner[:,:,i] + i * (order+1)^2;
        ind_inner[:,:,i] += i * (order + 1)^2;
    end
end

e = 1;
idx =  Mesh.get_node_indices(self, e, order);
ind1 = repeat(idx, NP, 1);
ind2 = reshape(repeat(idx', NP, 1), NPNP, 1);
st = (e - 1) * NPNP + 1;
en = e * NPNP;

I[st:en] = ind1;
J[st:en] = ind2;
pts =  Mesh.element_nodes(self, e, refel);
(detJac, Jac) = Mesh.geometric_factors(self, refel, pts);

eMat = Mesh.element_mass(self, e, refel, detJac);
# function element_mass(self, eid, refel, J)
#   # element mass matrix
#   Md = refel.W .* J ;
#   Mds = Md[:,1]
#   Me = refel.Q' * diagm(Mds) * refel.Q;
#   return Me
# end
mass_val[st:en] = eMat[:];

eMat = Mesh.element_stiffness(self, e, refel, detJac, Jac);
stiff_val[st:en] = eMat[:];

eMat_inner_inv = inv(eMat[ind_inner[:],ind_inner[:]]);
eMat_inv = diagm(diag(eMat, 0));

eMat_inv[ind_inner[:],ind_inner[:]] =  eMat_inner_inv;
inv_stiff_val[st:en] = eMat_inv[:];

# mats = square_stokes()
# domain = 1;
#
#   A = mats["A"]; B = mats["B"]; Bx = mats["Bx"]; By = mats["By"];
#   f = mats["f"]; g = mats["g"]; xy = mats["xy"]; xyp = mats["xyp"];
#   bound = mats["bound"]; x = mats["x"]; y = mats["y"];
#   Q = mats["Q"];
#   if domain == 3
#     P = mats["P"]
#   else
#     P = zeros(size(A))
#   end
#
#   # boundary conditions
#   println("imposing (enclosed flow) boundary conditions ...")
#   (Ast, Bst, fst, gst) = flowbc(A, B, f, g, xy, bound, domain)
#   np = length(gst)
#   rhs = vec([fst; gst])
#
#   ## compute solution
#   K = [Ast Bst'; Bst spzeros(np, np)]
#
# u = linspace(1, 62, 62)
#
#   nv = size(Ast, 1); np = size(Q, 1); nu = nv/2
#   Agal = Ast[1:nu, 1:nu]
#   (mgdata, smooth_data, sweeps, stype, npre, npost, nc) = mg_diff(x, y, Agal)
#
#   mparams = {
#     "nv" => nv,
#     "Q" => Q,
#     "mgdata" => mgdata,
#     "smooth_data" => smooth_data,
#     "nc" => nc,
#     "npre" => npre,
#     "npost" => npost,
#     "sweeps" => sweeps
#   }
#
#   # block GMG preconditioner
#   M = u -> m_st_mg(u, mparams)
#   w = M(u)

  # wa = Ast*u[1:50]
  # wb = Bst * u[1:50]
  # w2 = K*u
