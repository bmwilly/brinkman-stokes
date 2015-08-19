
module Mesh
include("Refel.jl")
include("Basis.jl")
include("Tensor.jl")
include("BaseCustom.jl")

import PyPlot
using Distances
export set_coeff

type Hexmesh
	dim#==2	TODO=#
	nelems#==[8 8]	TODO=#
	order
	Xf # transform
	# problem specific
	coeff

	function Hexmesh(nelems::Any, X::Any)
		mesh = new()
		mesh.dim = length(nelems)
		mesh.nelems = nelems
		mesh.Xf = X
		f(x,y,z=[])=1
		mesh.coeff = f
		mesh.order=[]
		return mesh
	end
end
type C
	num_nodes;
	num_elements;
	num_bdy_nodes
	num_int_elements;
	num_bdy_elements;
	nnz;
	function C()
		c = new()
		return c;
	end
end
type D
	rx;
	ry;
	rz;
	sx;
	sy;
	sz;
	tx;
	ty;
	tz;
	function D()
		d= new()
		return d
	end
end
function plot(self)

	# display the mesh. Needs X.
	# PyPlot.figure(1)
	def_color = (31/256,171/256,226/256);   # default color of grid
	lw = 1;                         # default line width of grid

	if (self.dim == 2 )
	  	(x,y) = ndgrid(0:1/self.nelems[1]:1.0, 0:1/self.nelems[2]:1.0);
	  	pts = [x[:] y[:]];
	  	coords = self.Xf(pts);
	  	PyPlot.plot(coords[:,1], coords[:,2], c="black", marker="o", linestyle="None", hold="True");
	    # hold on;
	    x = reshape(coords[:,1], self.nelems[1]+1, self.nelems[2]+1);
	    y = reshape(coords[:,2], self.nelems[1]+1, self.nelems[2]+1);
	    PyPlot.plot(x,y, c=def_color, linewidth=lw);
	    PyPlot.plot(x',y', c=def_color, linewidth=lw);
	    # boundary test
	    idx = get_boundary_node_indices(self, 1);
	    PyPlot.plot(coords[idx,1], coords[idx,2], c="r", marker="o", linestyle="None");
	    PyPlot.title(string("Quad Mesh - ", self.nelems[1], "x", self.nelems[2]))
	    # axis square
	else
	    # 2 xy planes
	    (x,y) = ndgrid( 0:1/self.nelems[1]:1.0, 0:1/self.nelems[2]:1.0 );

	    # z = 0
	    z = zeros (size(x));
	    pts = [x[:] y[:] z[:]];
	    coords = self.Xf(pts);
	    x0 = reshape(coords[:,1], size(x));
	    y0 = reshape(coords[:,2], size(x));
	    z0 = reshape(coords[:,3], size(x));
	    PyPlot.surf(x0, y0, z0, alpha=0.5);

	    # hold on;
	    x1 = reshape(coords[:,1], self.nelems[1]+1, self.nelems[2]+1);
	    y1 = reshape(coords[:,2], self.nelems[1]+1, self.nelems[2]+1);
	    z1 = reshape(coords[:,3], self.nelems[1]+1, self.nelems[2]+1);

	    PyPlot.plot3D(x1[:],y1[:],z1[:], c="k", linestyle="-");
	    PyPlot.plot3D(x1'[:],y1'[:],z1'[:], c="k",  linestyle="-");
	    # z = 1
	    z = ones  (size(x));
	    pts = [x[:] y[:] z[:]];
	    coords = self.Xf(pts);
	    PyPlot.surf(reshape(coords[:,1], size(x)), reshape(coords[:,2], size(x)), reshape(coords[:,3], size(x)), alpha=0.5);
	    # hold on;
	    x1 = reshape(coords[:,1], self.nelems[1]+1, self.nelems[2]+1);
	    y1 = reshape(coords[:,2], self.nelems[1]+1, self.nelems[2]+1);
	    z1 = reshape(coords[:,3], self.nelems[1]+1, self.nelems[2]+1);
	    PyPlot.plot3D(x1[:],y1[:],z1[:], c="k", linestyle="-");
	    PyPlot.plot3D(x1'[:],y1'[:],z1'[:], c="k", linestyle="-");
	    #--------------------------------------------------------------------
	    # 2 yz planes
	    (y,z) = ndgrid( 0:1/self.nelems[2]:1.0, 0:1/self.nelems[3]:1.0 );
	    # x = 0
	    x = zeros (size(y));
	    pts = [x[:] y[:] z[:]];
	    coords = self.Xf(pts);
	    PyPlot.surf(reshape(coords[:,1], size(x)), reshape(coords[:,2], size(x)), reshape(coords[:,3], size(x)), alpha=0.5);
	    # hold on;
	    x1 = reshape(coords[:,1], self.nelems[2]+1, self.nelems[3]+1);
	    y1 = reshape(coords[:,2], self.nelems[2]+1, self.nelems[3]+1);
	    z1 = reshape(coords[:,3], self.nelems[2]+1, self.nelems[3]+1);
	    PyPlot.plot3D(x1[:],y1[:],z1[:], c="k", linestyle="-");
	    PyPlot.plot3D(x1'[:],y1'[:],z1'[:], c="k",  linestyle="-");
	    # x = 1
	    x = ones (size(y));
	    pts = [x[:] y[:] z[:]];
	    coords = self.Xf(pts);
	    PyPlot.surf(reshape(coords[:,1], size(x)), reshape(coords[:,2], size(x)), reshape(coords[:,3], size(x)), alpha=0.5);
	    # hold on;
	    x1 = reshape(coords[:,1], self.nelems[2]+1, self.nelems[3]+1);
	    y1 = reshape(coords[:,2], self.nelems[2]+1, self.nelems[3]+1);
	    z1 = reshape(coords[:,3], self.nelems[2]+1, self.nelems[3]+1);
	    PyPlot.plot3D(x1[:],y1[:],z1[:], c="k", linestyle="-");
	    PyPlot.plot3D(x1'[:],y1'[:],z1'[:], c="k",  linestyle="-");
	    #--------------------------------------------------------------------
	    # 2 xz planes
	    (x,z) = ndgrid( 0:1/self.nelems[1]:1.0, 0:1/self.nelems[3]:1.0 );
	    # y = 0
	    y = zeros (size(x));
	    pts = [x[:] y[:] z[:]];
	    coords = self.Xf(pts);
	    PyPlot.surf(reshape(coords[:,1], size(x)), reshape(coords[:,2], size(x)), reshape(coords[:,3], size(x)), alpha=0.5);
	    # hold on;
	    x1 = reshape(coords[:,1], self.nelems[1]+1, self.nelems[3]+1);
	    y1 = reshape(coords[:,2], self.nelems[1]+1, self.nelems[3]+1);
	    z1 = reshape(coords[:,3], self.nelems[1]+1, self.nelems[3]+1);
	    PyPlot.plot3D(x1[:],y1[:],z1[:], c="k", linestyle="-");
	    PyPlot.plot3D(x1'[:],y1'[:],z1'[:], c="k",  linestyle="-");
	    # y = 1
	    y = ones (size(x));
	    pts = [x[:] y[:] z[:]];
	    coords = self.Xf(pts);
	    PyPlot.surf(reshape(coords[:,1], size(x)), reshape(coords[:,2], size(x)), reshape(coords[:,3], size(x)), alpha=0.5);
	    # hold on;
	    x1 = reshape(coords[:,1], self.nelems[1]+1, self.nelems[3]+1);
	    y1 = reshape(coords[:,2], self.nelems[1]+1, self.nelems[3]+1);
	    z1 = reshape(coords[:,3], self.nelems[1]+1, self.nelems[3]+1);
	    PyPlot.plot3D(x1[:],y1[:],z1[:], c="k", linestyle="-");
	    PyPlot.plot3D(x1'[:],y1'[:],z1'[:], c="k",  linestyle="-");

	    # pretty views etc
	    # view(3); axis equal;
	    # title(['Hex Mesh ' num2str(self.nelems[1]) 'x' num2str(self.nelems[2]) 'x' num2str(self.nelems[3])])
	end
	# set (gcf, 'renderer', 'opengl');
	# cameratoolbar('show');
	# cameratoolbar('setmode', 'orbit');
end
function plot_fx(self, fx)
# display the mesh. Needs X.
# hFig = figure(1);
# set(gcf,'PaperPositionMode','auto')
# set(hFig, 'Position', [200 200 800 800])

if (self.dim == 2 )
    scl=8;
    (x,y) = ndgrid(0:1/(scl*self.nelems[1]):1.0, 0:1/(scl*self.nelems[2]):1.0);
    pts = [x[:] y[:]];
    coords = self.Xf(pts);

    ci = arrayfun( fx, coords[:,1], coords[:,2] );
    # surf(reshape(coords[:,1], size(x)), ...
    #    reshape(coords[:,2], size(x)), ...
    #    zeros(size(x)), reshape(ci, size(x)), ...
    #    'EdgeColor','none','LineStyle','none' ...
    #    );
	PyPlot.surf(reshape(coords[:,1], size(x)), reshape(coords[:,2], size(x)), reshape(coords[:,3], size(x)), alpha=0.5);

    # hold on;
    (x,y) = ndgrid(0:1/self.nelems[1]:1.0, 0:1/self.nelems[2]:1.0);
    pts = [x[:] y[:]];
    coords = self.Xf(pts);
    x = reshape(coords[:,1], self.nelems[1]+1, self.nelems[2]+1);
    y = reshape(coords[:,2], self.nelems[1]+1, self.nelems[2]+1);
    PyPlot.plot(x,y, c="k", linestyle="-");
    PyPlot.plot(x',y', c="k", linestyle="-");

    # axis square;
    # view(0,90);
    # colorbar;
else
    # will draw 6 planes ...
    #--------------------------------------------------------------------
    # 2 xy planes
    scl=8;
    (x,y) = ndgrid( 0:1/(scl*self.nelems[1]):1.0, 0:1/(scl*self.nelems[2]):1.0 );
    # z = 0
    z = zeros (size(x));
    pts = [x[:] y[:] z[:]];
    coords = self.Xf(pts);
    ci = arrayfun( fx, coords[:,1], coords[:,2], coords[:,3] );
    # surf(reshape(coords[:,1], size(x)), ...
    #    reshape(coords[:,2], size(x)), ...
    #    reshape(coords[:,3], size(x)), ...
    #    reshape(ci, size(x)), ...
    #    'EdgeColor','none','LineStyle','none' ...
    #    );
	PyPlot.surf(reshape(coords[:,1], size(x)), reshape(coords[:,2], size(x)), reshape(coords[:,3], size(x)), alpha=0.5);

    # hold on;
    (x,y) = ndgrid( 0:1/self.nelems[1]:1.0, 0:1/self.nelems[2]:1.0 );
    z = zeros (size(x));
    pts = [x[:] y[:] z[:]];
    coords = self.Xf(pts);
    x1 = reshape(coords[:,1], self.nelems[1]+1, self.nelems[2]+1);
    y1 = reshape(coords[:,2], self.nelems[1]+1, self.nelems[2]+1);
    z1 = reshape(coords[:,3], self.nelems[1]+1, self.nelems[2]+1);
    PyPlot.plot3D(x1[:],y1[:],z1[:], c="k", linestyle="-");
    PyPlot.plot3D(x1'[:],y1'[:],z1'[:], c="k", linestyle="-");
    # z = 1
    (x,y) = ndgrid( 0:1/(scl*self.nelems[1]):1.0, 0:1/(scl*self.nelems[2]):1.0 );
    z = ones  (size(x));
    pts = [x[:] y[:] z[:]];
    coords = self.Xf(pts);
    ci = arrayfun( fx, coords[:,1], coords[:,2], coords[:,3] );
    PyPlot.surf(reshape(coords[:,1], size(x)), reshape(coords[:,2], size(x)), reshape(coords[:,3], size(x)), alpha=0.5);

    # surf(reshape(coords[:,1], size(x)), ...
    #    reshape(coords[:,2], size(x)), ...
    #    reshape(coords[:,3], size(x)), ...
    #    reshape(ci, size(x)), ...
    #    'EdgeColor','none','LineStyle','none' ...
    #     ); # 'FaceColor', 'interp', 'FaceLighting', 'phong'
    # hold on;
    (x,y) = ndgrid( 0:1/self.nelems[1]:1.0, 0:1/self.nelems[2]:1.0 );
    z = ones  (size(x));
    pts = [x[:] y[:] z[:]];
    coords = self.Xf(pts);
    x1 = reshape(coords[:,1], self.nelems[1]+1, self.nelems[2]+1);
    y1 = reshape(coords[:,2], self.nelems[1]+1, self.nelems[2]+1);
    z1 = reshape(coords[:,3], self.nelems[1]+1, self.nelems[2]+1);
    PyPlot.plot3D(x1[:],y1[:],z1[:], c="k", linestyle="-");
    PyPlot.plot3D(x1'[:],y1'[:],z1'[:], c="k", linestyle="-");
    #--------------------------------------------------------------------
    # 2 yz planes
    (y,z) = ndgrid( 0:1/(scl*self.nelems[2]):1.0, 0:1/(scl*self.nelems[3]):1.0 );
    # x = 0
    x = zeros (size(y));
    pts = [x[:] y[:] z[:]];
    coords = self.Xf(pts);
    ci = arrayfun( fx, coords[:,1], coords[:,2], coords[:,3] );
    # surf(reshape(coords[:,1], size(x)), ...
    #    reshape(coords[:,2], size(x)), ...
    #    reshape(coords[:,3], size(x)), ...
    #    reshape(ci, size(x)), ...
    #    'EdgeColor','none','LineStyle','none' ...
    #    );
    PyPlot.surf(reshape(coords[:,1], size(x)), reshape(coords[:,2], size(x)), reshape(coords[:,3], size(x)), alpha=0.5);

    # hold on;
    (y,z) = ndgrid( 0:1/self.nelems[2]:1.0, 0:1/self.nelems[3]:1.0 );
    # x = 0
    x = zeros (size(y));
    pts = [x[:] y[:] z[:]];
    coords = self.Xf(pts);
    x1 = reshape(coords[:,1], self.nelems[2]+1, self.nelems[3]+1);
    y1 = reshape(coords[:,2], self.nelems[2]+1, self.nelems[3]+1);
    z1 = reshape(coords[:,3], self.nelems[2]+1, self.nelems[3]+1);
    PyPlot.plot3D(x1[:],y1[:],z1[:], c="k", linestyle="-");
    PyPlot.plot3D(x1'[:],y1'[:],z1'[:], c="k", linestyle="-");
    x = 1
    (y,z) = ndgrid( 0:1/(scl*self.nelems[2]):1.0, 0:1/(scl*self.nelems[3]):1.0 );
    x = ones (size(y));
    pts = [x[:] y[:] z[:]];
    coords = self.Xf(pts);
    ci = arrayfun( fx, coords[:,1], coords[:,2], coords[:,3] );
    # surf(reshape(coords[:,1], size(x)), ...
    #    reshape(coords[:,2], size(x)), ...
    #    reshape(coords[:,3], size(x)), ...
    #    reshape(ci, size(x)), ...
    #    'EdgeColor','none','LineStyle','none' ...
    #    );
    PyPlot.surf(reshape(coords[:,1], size(x)), reshape(coords[:,2], size(x)), reshape(coords[:,3], size(x)), alpha=0.5);

    # hold on;
    (y,z) = ndgrid( 0:1/self.nelems[2]:1.0, 0:1/self.nelems[3]:1.0 );
    x = ones (size(y));
    pts = [x[:] y[:] z[:]];
    coords = self.Xf(pts);
    x1 = reshape(coords[:,1], self.nelems[2]+1, self.nelems[3]+1);
    y1 = reshape(coords[:,2], self.nelems[2]+1, self.nelems[3]+1);
    z1 = reshape(coords[:,3], self.nelems[2]+1, self.nelems[3]+1);
    PyPlot.plot3D(x1[:],y1[:],z1[:], c="k", linestyle="-");
    PyPlot.plot3D(x1'[:],y1'[:],z1'[:], c="k", linestyle="-");
    #--------------------------------------------------------------------
    # 2 xz planes
    (x,z) = ndgrid( 0:1/(scl*self.nelems[1]):1.0, 0:1/(scl*self.nelems[3]):1.0 );
    # y = 0
    y = zeros (size(x));
    pts = [x[:] y[:] z[:]];
    coords = self.Xf(pts);
    ci = arrayfun( fx, coords[:,1], coords[:,2], coords[:,3] );
    # surf(reshape(coords[:,1], size(x)), ...
    #    reshape(coords[:,2], size(x)), ...
    #    reshape(coords[:,3], size(x)), ...
    #    reshape(ci, size(x)) ...
    #    );
    PyPlot.surf(reshape(coords[:,1], size(x)), reshape(coords[:,2], size(x)), reshape(coords[:,3], size(x)), alpha=0.5);

    # hold on;
    (x,z) = ndgrid( 0:1/self.nelems[1]:1.0, 0:1/self.nelems[3]:1.0 );
    y = zeros (size(x));
    pts = [x[:] y[:] z[:]];
    coords = self.Xf(pts);
    x1 = reshape(coords[:,1], self.nelems[1]+1, self.nelems[3]+1);
    y1 = reshape(coords[:,2], self.nelems[1]+1, self.nelems[3]+1);
    z1 = reshape(coords[:,3], self.nelems[1]+1, self.nelems[3]+1);
    PyPlot.plot3D(x1[:],y1[:],z1[:], c="k", linestyle="-");
    PyPlot.plot3D(x1'[:],y1'[:],z1'[:], c="k", linestyle="-");
    # y = 1
    (x,z) = ndgrid( 0:1/(scl*self.nelems[1]):1.0, 0:1/(scl*self.nelems[3]):1.0 );
    y = ones (size(x));
    pts = [x[:] y[:] z[:]];
    coords = self.Xf(pts);
    ci = arrayfun( fx, coords[:,1], coords[:,2], coords[:,3] );
    # surf(reshape(coords[:,1], size(x)), ...
    #    reshape(coords[:,2], size(x)), ...
    #    reshape(coords[:,3], size(x)), ...
    #    reshape(ci, size(x)) ...
    #    );
    PyPlot.surf(reshape(coords[:,1], size(x)), reshape(coords[:,2], size(x)), reshape(coords[:,3], size(x)), alpha=0.5);

    # hold on;
    (x,z) = ndgrid( 0:1/self.nelems[1]:1.0, 0:1/self.nelems[3]:1.0 );
    y = ones (size(x));
    pts = [x[:] y[:] z[:]];
    coords = self.Xf(pts);
    x1 = reshape(coords[:,1], self.nelems[1]+1, self.nelems[3]+1);
    y1 = reshape(coords[:,2], self.nelems[1]+1, self.nelems[3]+1);
    z1 = reshape(coords[:,3], self.nelems[1]+1, self.nelems[3]+1);
    PyPlot.plot3D(x1[:],y1[:],z1[:], c="k", linestyle="-");
    PyPlot.plot3D(x1'[:],y1'[:],z1'[:], c="k", linestyle="-");

    # pretty views etc
    # view(3); axis square
    # view(150,40);
    # # colorbar;
end

# set (gcf, 'renderer', 'opengl');
# cameratoolbar('show');
# cameratoolbar('setmode', 'orbit');
end
	function set_order(self, order)
		if isempty(self.order)
			self.order = order;
		else
			assert (order == self.order);
		end
	end
	function set_coeff(self, coeff)
		if ( typeof(coeff) == ASCIIString )
			# is a string, so convert into a function
#			syms x y z;
#TODO			expr = ['matlabFunction(' coeff ')'];
			self.coeff = eval(expr);
		else
			self.coeff = coeff;
		end
	end
	function assemble_mass(self, order)
		set_order(self, order);
		# assemble the mass matrix
		refel = Refel( self.dim, order );
		dof = prod([self.nelems...]*order + 1);
		ne = prod([self.nelems...]);
		# storage for indices and values
		NP = (order+1)^self.dim;
		NPNP = NP * NP;
		eM = zeros(NP, NP);
		I = zeros(ne * NPNP, 1);
		J = zeros(ne * NPNP, 1);
		val = zeros(ne * NPNP, 1);
		# loop over elements
		for e=1:ne
			pts =  element_nodes(self, e, refel);
			(detJac,D) = geometric_factors(self, refel, pts);
			idx =  get_node_indices (self, e, order);
			eM = element_mass(self, e, refel, detJac);
			ind1 = repmat(idx,NP,1);
			ind2 = reshape(repmat(idx',NP,1),NPNP,1);
			st = (e-1)*NPNP+1;
			en = e*NPNP;
			I[st:en] = ind1;
			J[st:en] = ind2;
			val[st:en] = eM[:]
		end
		M = sparse(int64(I[:]),int64(J[:]),val[:],dof,dof);
		return M;
	end
	function assemble_stiffness(self, order)
		set_order(self, order);
		# assemble the stiffness matrix
		refel = Refel( self.dim, order );
		dof = prod([self.nelems...]*order + 1);
		ne = prod([self.nelems...]);
		# storage for indices and values
		NP = (order+1)^self.dim;
		NPNP = NP * NP;
		I = zeros(ne * NPNP, 1);
		J = zeros(ne * NPNP, 1);
		stiff_val = zeros(ne * NPNP, 1);
		# loop over elements
		for e=1:ne
			idx =  get_node_indices (self, e, order);
			ind1 = repmat(idx,NP,1);
			ind2 = reshape(repmat(idx',NP,1),NPNP,1);
			st = (e-1)*NPNP+1;
			en = e*NPNP;
			I[st:en] = ind1;
			J[st:en] = ind2;
			pts =  element_nodes(self, e, refel);
			(detJac, Jac) = geometric_factors(self, refel, pts);
			eMat = element_stiffness(self, e, refel, detJac, Jac);
			stiff_val[st:en] = eMat[:];
		end
		return K
	end
	function assemble_poisson(self, order)
		set_order(self,order);
		# assemble the mass matrix
		refel = Refel( self.dim, order );
		dof = prod([self.nelems...]*order + 1);
		ne = prod([self.nelems...]);
		# storage for indices and values
		NP = (order+1)^self.dim;
		NPNP = NP * NP;

		I = zeros(ne * NPNP, 1);
		J = zeros(ne * NPNP, 1);
		mass_val = zeros(ne * NPNP, 1);
		stiff_val = zeros(ne * NPNP, 1);
		inv_stiff_val = zeros(ne * NPNP, 1);
		ind_inner1D = repmat((2:order), 1, order-1);
		if self.dim == 2

			ind_inner = ind_inner1D + (order+1) * (ind_inner1D'-1);
		else
			ind_inner = ind_inner1D + (order+1) * (ind_inner1D'-1);
			# ind_inner = repmat(ind_inner, [1,1,order-1]);
			ind_inner = repeat(ind_inner, outer = [1,1,order-1]);
			for i = 1:order-1
				# ind_inner[:,:,i] = ind_inner[:,:,i] + i * (order+1)^2;
				ind_inner[:,:,i] += i * (order+1)^2;
			end
		end

		# loop over elements
		for e=1:ne
			idx =  get_node_indices(self, e, order);
			ind1 = repmat(idx,NP,1);
			ind2 = reshape(repmat(idx',NP,1),NPNP,1);
			st = (e-1)*NPNP+1;
			en = e*NPNP;

			I[st:en] = ind1;
			J[st:en] = ind2;
			pts =  element_nodes(self, e, refel);
			(detJac, Jac) = geometric_factors(self, refel, pts);

			eMat = element_mass(self, e, refel, detJac);
			mass_val[st:en] = eMat[:];

			eMat = element_stiffness(self, e, refel, detJac, Jac);
			stiff_val[st:en] = eMat[:];

			eMat_inner_inv = inv(eMat[ind_inner[:],ind_inner[:]]);
			eMat_inv = diagm(diag(eMat,0));

			eMat_inv[ind_inner[:],ind_inner[:]] =  eMat_inner_inv;
			inv_stiff_val[st:en] = eMat_inv[:];
		end

		Iv=int64(I[:]);
		Jv=int64(J[:]);
		mv=mass_val[:];

		M = sparse(Iv,Jv,mv,dof,dof);
		# zero dirichlet bdy conditions
		bdy = get_boundary_node_indices(self, order);
		ii = ismember(I,bdy);
		jj = ismember(J,bdy);

		stiff_val = stiff_val.*(int(!bool(ii))).*(int(!bool(jj)));
		inv_stiff_val = inv_stiff_val.*(int(!bool(ii))).*(int(!bool(jj)));

		I = [I; bdy];
		J = [J; bdy];
		stiff_val = [stiff_val; ones(length(bdy), 1)];
		inv_stiff_val = [inv_stiff_val; ones(length(bdy), 1)];
		Iv=int64(I[:]);
		Jv=int64(J[:]);
		sv=stiff_val[:];
		isv=inv_stiff_val[:];

		K = sparse(Iv,Jv,sv,dof,dof);
		iK = sparse(Iv,Jv,isv,dof,dof);
		ebdy = get_element_boundary_node_indices(self, order);
		iKebdry = diag(full(iK[ebdy,ebdy]),0)
		if countnz(iKebdry) > 0
      	iK[ebdy,ebdy] = diagm(1./iKebdry)
	  end
		return K, M, iK
	end
	function assemble_poisson_brinkman(self, order, centers)
		set_order(self,order);
		# assemble the mass matrix
		refel = Refel( self.dim, order );
		dof = prod([self.nelems...]*order + 1);
		ne = prod([self.nelems...]);
		# storage for indices and values
		NP = (order+1)^self.dim;
		NPNP = NP * NP;

		I = zeros(ne * NPNP, 1);
		J = zeros(ne * NPNP, 1);
		mass_val = zeros(ne * NPNP, 1);
		stiff_val = zeros(ne * NPNP, 1);
		perm_val = zeros(dof, 1);
		# inv_stiff_val = zeros(ne * NPNP, 1);
		# ind_inner1D = repmat((2:order), 1, order-1);
		# if self.dim == 2
		#
		# 	ind_inner = ind_inner1D + (order+1) * (ind_inner1D'-1);
		# else
		# 	ind_inner = ind_inner1D + (order+1) * (ind_inner1D'-1);
		# 	# ind_inner = repmat(ind_inner, [1,1,order-1]);
		# 	ind_inner = repeat(ind_inner, outer = [1,1,order-1]);
		# 	for i = 1:order-1
		# 		# ind_inner[:,:,i] = ind_inner[:,:,i] + i * (order+1)^2;
		# 		ind_inner[:,:,i] += i * (order+1)^2;
		# 	end
		# end

		# loop over elements
		for e=1:ne
			idx =  get_node_indices(self, e, order);
			ind1 = repmat(idx,NP,1);
			ind2 = reshape(repmat(idx',NP,1),NPNP,1);
			st = (e-1)*NPNP+1;
			en = e*NPNP;

			I[st:en] = ind1;
			J[st:en] = ind2;
			pts =  element_nodes(self, e, refel);
			(detJac, Jac) = geometric_factors(self, refel, pts);
			brinkman_pts = brinkman_tensor(pts, centers);

			eMat = element_mass_brinkman(self, e, refel, detJac, brinkman_pts);
			mass_val[st:en] = eMat[:];

			eMat = element_stiffness_brinkman(self, e, refel, detJac, Jac, brinkman_pts);
			stiff_val[st:en] = eMat[:];

			perm_val[idx] = brinkman_pts[:];

			# eMat_inner_inv = inv(eMat[ind_inner[:],ind_inner[:]]);
			# eMat_inv = diagm(diag(eMat,0));
			#
			# eMat_inv[ind_inner[:],ind_inner[:]] =  eMat_inner_inv;
			# inv_stiff_val[st:en] = eMat_inv[:];
		end

		Iv=int64(I[:]);
		Jv=int64(J[:]);
		mv=mass_val[:];

		M = sparse(Iv,Jv,mv,dof,dof);
		# zero dirichlet bdy conditions
		bdy = get_boundary_node_indices(self, order);
		ii = ismember(I,bdy);
		jj = ismember(J,bdy);

		stiff_val = stiff_val.*(int(!bool(ii))).*(int(!bool(jj)));
		# inv_stiff_val = inv_stiff_val.*(int(!bool(ii))).*(int(!bool(jj)));

		I = [I; bdy];
		J = [J; bdy];
		stiff_val = [stiff_val; ones(length(bdy), 1)];
		# inv_stiff_val = [inv_stiff_val; ones(length(bdy), 1)];
		Iv=int64(I[:]);
		Jv=int64(J[:]);
		sv=stiff_val[:];
		# isv=inv_stiff_val[:];

		K = sparse(Iv,Jv,sv,dof,dof);
		# iK = sparse(Iv,Jv,isv,dof,dof);
		# ebdy = get_element_boundary_node_indices(self, order);
		# iKebdry = diag(full(iK[ebdy,ebdy]),0)
		# if countnz(iKebdry) > 0
    #   	iK[ebdy,ebdy] = diagm(1./iKebdry)
	  # end
		# return K, M, iK

		kappa = perm_val[:];

		return K,M,kappa
	end
	function assemble_permeability(self, order, centers)
		set_order(self, order);
		refel = Refel(self.dim, order);
		dof = prod([self.nelems...]*order + 1);
		ne = prod([self.nelems...]);
		perm_val = zeros(dof, 1);
		# loop over elements
		for e = 1:ne
			idx = get_node_indices(self, e, order);
			pts = element_nodes(self, e, refel);
			brinkman_pts = brinkman_tensor(pts, centers);
			perm_val[idx] = brinkman_pts[:];
		end
		kappa = perm_val[:];
	end
	function assemble_rhs(self, fx, order)
		set_order(self, order)
		refel = Refel(self.dim, order)
		dof = prod([self.nelems...]*order + 1)
		ne = prod([self.nelems...])
		f = zeros(dof,1)
		# loop over elements
		for e=1:ne
			idx =  get_node_indices (self, e, order)
			pts =  element_nodes(self, e, refel)
			(J,D) = geometric_factors(self, refel, pts)
			gpts = element_gauss(self, e, refel)
			if (self.dim == 2)
				fd = fx(gpts[:,1], gpts[:,2] )
			else
				fd = fx(gpts[:,1], gpts[:,2], gpts[:,3] )
			end
			Jd = refel.W .* J .* fd
			f[idx] = f[idx] + refel.Q' * Jd
		end
		return f
	end
	function assemble_interpolation(self, order)
		# assemble prolongation operator from coarse (self) to fine mesh
		refel = Refel( self.dim, self.order );

		if ( order == self.order )
			dof_coarse = prod([self.nelems...] * self.order + 1);
			dof_fine   = prod(2*[self.nelems...] * self.order + 1);
			NP_c = (self.order+1)^self.dim;
			NP_f = (2*self.order+1)^self.dim;
			Pe = refel.Ph;
		else
			assert (order == 2*self.order);
			NP_c = (self.order+1)^self.dim;
			NP_f = (order+1)^self.dim;
			dof_coarse = prod([self.nelems...] * self.order + 1);
			dof_fine   = prod([self.nelems...] * order + 1);
			Pe = refel.Pp;
		end
		ne  = prod([self.nelems...]);
		# storage for indices and values
		NPNP = NP_c * NP_f;
		I = zeros(ne * NPNP, 1);
		J = zeros(ne * NPNP, 1);
		val = zeros(ne * NPNP, 1);

		for e=1:ne
			(idx_c, idx_f) = get_interpolation_indices (self,e);
			ind1 = repmat(idx_f,NP_c,1);
			ind2 = reshape(repmat(idx_c',NP_f,1),NPNP,1);
			st = (e-1)*NPNP+1;
			en = e*NPNP;
			I[st:en] = ind1;
			J[st:en] = ind2;

			val[st:en] = reshape(Pe,size(Pe)[1].*size(Pe)[2],1)
		end
		u_ij = unique([I J],1);
		q=indunique([I J],1);
		u_val   = val[q];
		I = int64(u_ij[:,1]);
		J = int64(u_ij[:,2]);
		P = sparse(I[:],J[:],u_val[:],dof_fine,dof_coarse);
		return P;
	end
	function get_node_indices ( self, eid, order )
		# determine global node indices for a given element

		if ( self.dim == 2)
			(i,j) = ind2sub (self.nelems, eid);
			i_low   = (i-1)*order + 1;   i_high =  i*order + 1;
			j_low   = (j-1)*order + 1;   j_high =  j*order + 1;

			(i,j) = ndgrid(i_low:i_high, j_low:j_high);
			m=i[:]
			n=j[:]
			x=[1:length(m)]
			idx = sub2ind([self.nelems...]'*order + 1,m[x],n[x])
		else
			(i,j,k) = ind2sub (self.nelems, eid);

			i_low   = (i-1)*order + 1;   i_high =  i*order + 1;
			j_low   = (j-1)*order + 1;   j_high =  j*order + 1;
			k_low   = (k-1)*order + 1;   k_high =  k*order + 1;

			(i,j,k) = ndgrid(i_low:i_high, j_low:j_high, k_low:k_high);

			m=i[:]
			n=j[:]
			o=k[:]
			x=[1:length(m)]
			idx = sub2ind([self.nelems...]'*order + 1,m[x],n[x],o[x])
		end
		return idx
	end
	function get_linear_node_indices ( self, eid, order )
		# determine global node indices for a given element
		if ( self.dim == 2)
			(i,j) = ind2sub (self.nelems*order, eid);

			(i,j) = ndgrid(i:i+1, j:j+1);

			idx     = sub2ind ([self.nelems...]'*order + 1, i[:], j[:]);
		else
			(i,j,k) = ind2sub (self.nelems*order, eid);

			(i,j,k) = ndgrid(i:i+1, j:j+1, k:k+1);

			idx     = sub2ind ([self.nelems...]'*order + 1, i[:], j[:], k[:] );
		end
		return idx
	end
	function get_interpolation_indices ( self, eid )
		# determine global node indices for a given element
		if ( self.dim == 2)
			(i,j) = ind2sub (self.nelems, eid);

			i_low       = (i-1)*self.order + 1;   i_high =  i*self.order + 1;
			j_low       = (j-1)*self.order + 1;   j_high =  j*self.order + 1;
			(i,j)       = ndgrid(i_low:i_high, j_low:j_high);
			idx_coarse  = sub2ind ([self.nelems...]*self.order + 1, i[:], j[:]);

			(i,j)       = ndgrid(2*i_low-1:2*i_high-1, 2*j_low-1:2*j_high-1);
			idx_fine    = sub2ind (2*[self.nelems...]*self.order + 1, i[:], j[:]);
		else
			(i,j,k) = ind2sub (self.nelems, eid);

			i_low       = (i-1)*self.order + 1;   i_high =  i*self.order + 1;
			j_low       = (j-1)*self.order + 1;   j_high =  j*self.order + 1;
			k_low       = (k-1)*self.order + 1;   k_high =  k*self.order + 1;
			(i,j,k)     = ndgrid(i_low:i_high, j_low:j_high, k_low:k_high);
			idx_coarse  = sub2ind ([self.nelems...]*self.order + 1, i[:], j[:], k[:] );

			(i,j,k)     = ndgrid(2*i_low-1:2*i_high-1, 2*j_low-1:2*j_high-1, 2*k_low-1:2*k_high-1);
			idx_fine    = sub2ind (2*[self.nelems...]*self.order + 1, i[:], j[:], k[:] );
		end
		return idx_coarse, idx_fine
	end
	function get_boundary_node_indices(self, order)
		# function idx = get_boundary_node_indices(self, order)
		#    returns indices of boundary nodes, for setting
		#    boundary conditions
		if (self.dim == 2)
			(x,y) = ndgrid(1:self.nelems[1]*order+1,1:self.nelems[2]*order+1);

			idx = [ findin(x,1);findin(x,(self.nelems[1]*order+1));findin(y,1);findin(y,(self.nelems[2]*order+1))];
			idx = unique(sort(idx));
		else
			(x,y,z) = ndgrid(1:self.nelems[1]*order+1,1:self.nelems[2]*order+1,1:self.nelems[3]*order+1);

			idx = [ findin(x,1); findin(x,(self.nelems[1]*order+1)); findin(y,1); findin(y,(self.nelems[2]*order+1)); findin(z,1); findin(z,(self.nelems[3]*order+1))];

			idx = unique(sort(idx));
		end
		return idx
	end
	function get_element_boundary_node_indices(self, order)
		# function idx = get_element_boundary_node_indices(self, order)
		#    returns indices of element boundary nodes, for block
		#    Jacobi smoother
		if (self.dim == 2)
			(x,y) = ndgrid(1:self.nelems[1]*order+1,1:self.nelems[2]*order+1);
			idx = [ findin(mod(x,order),1);findin(mod(y,order),1);];
			idx = unique(sort(idx));
		else
			(x,y,z) = ndgrid(1:self.nelems[1]*order+1,1:self.nelems[2]*order+1,1:self.nelems[3]*order+1);

			idx = [ findin(mod(x,order),1); findin(mod(y,order),1); findin(mod(z,order),1);];

			idx = unique(sort(idx));
		end
		return idx
	end
	function element_mass(self, eid, refel, J)
		# element mass matrix
		Md = refel.W .* J ;
		Mds = Md[:,1]
		Me = refel.Q' * diagm(Mds) * refel.Q;
		return Me
	end
	function element_mass_brinkman(self, eid, refel, J, brinkman_pts)
		# element mass matrix for brinkman
		Md = refel.W .* J .* brinkman_pts;
		Mds = Md[:,1]
		Me = refel.Q' * diagm(Mds) * refel.Q;
		return Me
	end
	function brinkman_tensor(pts, centers)
		npts = length(pts[:,1]); nc = length(centers[:,1])
		brinkman_pts = ones(npts)

		# get euclidean distances between nodal points and centers of brinkman obstacles
		R = pairwise(Euclidean(), pts', centers')
		R = minimum(R, 2)
		# brinkman_pts[find(R .< 0.025)] = 1e6
		brinkman_pts[find(R .< 0.075)] = 1e6
		# brinkman_pts[find(R .< 0.2)] = 1e6
		# brinkman_pts = centers[1:npts]
		brinkman_pts
	end


	function element_stiffness(self, eid, r, J, D)
		# element mass matrix

		#             | Qx Qy Qz || rx ry rz |     | rx sx tx || Qx |
		#    Ke =                 | sx sy sz | J W | ry sy ty || Qy |
		#                         | tx ty tz |     | rz sz tz || Qz |
		gpts = element_gauss(self, eid, r);
		nn = length(J);

		factor = zeros(nn, 6);

		#             1  4  5
		# factor      4  2  6
		#             5  6  3

		if (self.dim == 2 )
			mu = map(self.coeff, gpts[:,1], gpts[:,2] );
			factor [:,1] = (D.rx.*D.rx + D.ry.*D.ry ) .* J .* r.W .* mu ; # d2u/dx^2
			factor [:,2] = (D.sx.*D.sx + D.sy.*D.sy ) .* J .* r.W .* mu ; # d2u/dy^2
			factor [:,3] = (D.rx.*D.sx + D.ry.*D.sy ) .* J .* r.W .* mu ; # d2u/dxdy

			Ke =   r.Qx' * diagm(factor[:,1]) * r.Qx + r.Qy' * diagm(factor[:,2]) * r.Qy + r.Qx' * diagm(factor[:,3]) * r.Qy + r.Qy' * diagm(factor[:,3]) * r.Qx ;
		else
			mu = map(self.coeff, gpts[:,1], gpts[:,2], gpts[:,3] );

			# first compute dj.w.J.J'
			factor [:,1] = (D.rx.*D.rx + D.ry.*D.ry + D.rz.*D.rz ) .* J .* r.W .* mu ; # d2u/dx^2
			factor [:,2] = (D.sx.*D.sx + D.sy.*D.sy + D.sz.*D.sz ) .* J .* r.W .* mu ; # d2u/dy^2
			factor [:,3] = (D.tx.*D.tx + D.ty.*D.ty + D.tz.*D.tz ) .* J .* r.W .* mu ; # d2u/dz^2

			factor [:,4] = (D.rx.*D.sx + D.ry.*D.sy + D.rz.*D.sz ) .* J .* r.W .* mu ; # d2u/dxdy
			factor [:,5] = (D.rx.*D.tx + D.ry.*D.ty + D.rz.*D.tz ) .* J .* r.W .* mu ; # d2u/dxdz
			factor [:,6] = (D.sx.*D.tx + D.sy.*D.ty + D.sz.*D.tz ) .* J .* r.W .* mu ; # d2u/dydz

			Ke =   r.Qx' * diagm(factor[:,1]) * r.Qx + r.Qy' * diagm(factor[:,2]) * r.Qy + r.Qz' * diagm(factor[:,3]) * r.Qz + r.Qx' * diagm(factor[:,4]) * r.Qy + r.Qy' * diagm(factor[:,4]) * r.Qx + r.Qx' * diagm(factor[:,5]) * r.Qz + r.Qz' * diagm(factor[:,5]) * r.Qx + r.Qz' * diagm(factor[:,6]) * r.Qy + r.Qy' * diagm(factor[:,6]) * r.Qz ;
		end
		return Ke
	end
	function element_stiffness_brinkman(self, eid, r, J, D, brinkman_pts)
		eMass = element_mass_brinkman(self, eid, r, J, brinkman_pts)
		eStiff = element_stiffness(self, eid, r, J, D)
		return eStiff + eMass
	end

	function geometric_factors( self, refel, pts )
		# change to using Qx etc ?
		if (refel.dim == 1)
			xr  = refel.Dg * pts;
			J = xr;
		elseif (refel.dim == 2)
			(xr, xs) = Tensor.grad2 (refel.Dg, pts[:,1]);
			(yr, ys) = Tensor.grad2 (refel.Dg, pts[:,2]);
			J = -xs.*yr + xr.*ys;
		else
			(xr, xs, xt) = Tensor.grad3 (refel.Dg, pts[:,1]);
			(yr, ys, yt) = Tensor.grad3 (refel.Dg, pts[:,2]);
			(zr, zs, zt) = Tensor.grad3 (refel.Dg, pts[:,3]);

			J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
		end
		D=Mesh.D();
		if (refel.dim == 1)
			D.rx = 1./J;
		elseif (refel.dim == 2)
			D.rx =  ys./J;
			D.sx = -yr./J;
			D.ry = -xs./J;
			D.sy =  xr./J;
		else
			D.rx =  (ys.*zt - zs.*yt)./J;
			D.ry = -(xs.*zt - zs.*xt)./J;
			D.rz =  (xs.*yt - ys.*xt)./J;

			D.sx = -(yr.*zt - zr.*yt)./J;
			D.sy =  (xr.*zt - zr.*xt)./J;
			D.sz = -(xr.*yt - yr.*xt)./J;

			D.tx =  (yr.*zs - zr.*ys)./J;
			D.ty = -(xr.*zs - zr.*xs)./J;
			D.tz =  (xr.*ys - yr.*xs)./J;
		end
		return J,D
	end
	function linear_element_nodes(self, elem, order)

		if (self.dim == 2)
			(i,j) = ind2sub (self.nelems*order, elem);

			x1d = getGLLcoords(order, self.nelems(1));
			y1d = getGLLcoords(order, self.nelems(2));

			(x, y) = ndgrid(x1d(i:i+1), y1d(j:j+1));
			pts = [x[:] y[:]];
		else
			(i,j,k) = ind2sub (self.nelems*order, elem);

			x1d = getGLLcoords(order, self.nelems(1));
			y1d = getGLLcoords(order, self.nelems(2));
			z1d = getGLLcoords(order, self.nelems(3));

			(x, y, z) = ndgrid(x1d(i:i+1), y1d(j:j+1), z1d(k:k+1));
			pts = [x[:] y[:] z[:]];
		end

		coords = self.Xf(pts);
	end

	function element_nodes(self, elem, refel)
		h = 1./[self.nelems...]';
		if ( self.dim == 2)
			(i,j) = ind2sub (self.nelems, elem);
			idx = [i j];
		else
			(i,j,k) = ind2sub (self.nelems, elem);
			idx = [i j k];
		end
		p_mid = (idx - 0.5) .* h;
		p_gll = refel.r * 0.5 * h;
		nodes = p_mid .+ p_gll;
		if ( self.dim == 2)
			(x, y) = ndgrid(nodes[:,1], nodes[:,2]);
			pts = [x[:] y[:]];
		else
			(x, y, z) = ndgrid(nodes[:,1], nodes[:,2], nodes[:,3]);
			pts = [x[:] y[:] z[:]];
		end
		coords = self.Xf(pts);
	end

	function element_gauss(self, elem, refel)
		# function pts = element_gauss(self, elem, refel)
		# returns location of gauss coordinates of order
		# for element
		if (self.order == refel.N)
			h = 1./[self.nelems...]';

			if ( self.dim == 2)
				(i,j) = ind2sub (self.nelems, elem);
				idx = [i j];
			else
				(i,j,k) = ind2sub (self.nelems, elem);
				idx = [i j k];
			end
			p_mid = (idx - 0.5) .* h;
			p_gau = refel.g * 0.5 * h;
			nodes = p_mid .+ p_gau;
			if ( self.dim == 2)
				(x, y) = ndgrid(nodes[:,1], nodes[:,2]);
				pts = [x[:] y[:]];
			else
				(x, y, z) = ndgrid(nodes[:,1], nodes[:,2], nodes[:,3]);
				pts = [x[:] y[:] z[:]];
			end
		else
			assert(refel.N == 1);
			# ... get gll points ...
			if (self.dim == 2)
				(i,j) = ind2sub (tuple([self.nelems...]'*self.order), elem);

				x1d = getGLLcoords(self.order, self.nelems[1]);
				y1d = getGLLcoords(self.order, self.nelems[2]);

				xg = x1d(i) + (x1d(i+1) - x1d(i))*(refel.g + 1)*0.5;
				yg = y1d(j) + (y1d(j+1) - y1d(j))*(refel.g + 1)*0.5;

				(x, y) = ndgrid(xg, yg);

				pts = [x[:] y[:]];
			else
				(i,j,k) = ind2sub (tuple([self.nelems...]'*self.order), elem);

				x1d = getGLLcoords(self.order, self.nelems[1]);
				y1d = getGLLcoords(self.order, self.nelems[2]);
				z1d = getGLLcoords(self.order, self.nelems[3]);

				xg = x1d(i) + (x1d(i+1) - x1d(i))*(refel.g + 1)*0.5;
				yg = y1d(j) + (y1d(j+1) - y1d(j))*(refel.g + 1)*0.5;
				zg = z1d(k) + (z1d(k+1) - z1d(k))*(refel.g + 1)*0.5;

				(x, y, z) = ndgrid(xg, yg, zg);

				pts = [x[:] y[:] z[:]];
			end
		end
		coords = self.Xf(pts);
	end
	function assemble_poisson_linearized (self, order)
		set_order(self, order);

		refel = Refel ( self.dim, 1 );

		dof = prod ( [self.nelems...]*order + 1);
		ne  = prod ( [self.nelems...]*order ) ;

		# storage for indices and values
		NP = (1+1)^self.dim; # linear elements
		NPNP = NP * NP;
		# eMat = zeros(NP, NP);

		I = zeros(ne * NPNP, 1);
		J = zeros(ne * NPNP, 1);
		mass_val = zeros(ne * NPNP, 1);
		stiff_val = zeros(ne * NPNP, 1);

		# loop over elements
		for e=1:ne
			idx = get_linear_node_indices (self, e, order);

			ind1 = repmat(idx,NP,1);
			ind2 = reshape(repmat(idx',NP,1),NPNP,1);
			st = (e-1)*NPNP+1;
			en = e*NPNP;
			I[st:en] = ind1;
			J[st:en] = ind2;

			pts = linear_element_nodes(self, e, order);

			(detJac, Jac) = geometric_factors(self, refel, pts);

			eMat = element_mass(self, e, refel, detJac);
			mass_val[st:en] = eMat[:];

			eMat = element_stiffness(self, e, refel, detJac, Jac);
			stiff_val[st:en] = eMat[:];
		end
		M = sparse(int64(I[:]),int64(J[:]),mass_val[:],dof,dof);
		# zero dirichlet bdy conditions
		bdy = get_boundary_node_indices(self, order);

		ii = ismember(I,bdy);
		jj = ismember(J,bdy);

		stiff_val = stiff_val.*(int(!bool(ii))).*(int(!bool(jj)));
		I = [I; bdy];
		J = [J; bdy];
		stiff_val = [stiff_val; ones(length(bdy), 1)];

		K = sparse(int64(I[:]),int64(J[:]),stiff_val[:],dof,dof);
		return K, M
	end
	function getGLLcoords(order, elems)
		# function coords=getGLLcoords(order, elems)
		# returns location of gll coordinates of order
		# for elements in [0,1]

		fac = 1.0/(2*elems);

		# gll coordinates in [-1,1]
		x = Basis.gll (0,0,order)';

		x = (x + 1)*fac;


		coords = [];
		for i=1:elems
			y = x + (i-1)/elems;
			coords = [coords y[1:end-1]];
		end

		coords = [coords 1.0];
	end

	function getUniformCoords(order, elems)
		coords = linspace(0, 1, order*elems+1);
	end

	function getElementCenters(order, elems)
		# order is ignored ...
		nodes = linspace(0,1, elems+1);
		coords = 1/2*(nodes[1:end-1] + nodes[2:end]);
	end

	function stats(nelems, order)
		# function Ch = stats(nelems, order)
		#   given number of elements and the order,
		#   this function calculates different node
		#   stats for the mesh
		C=Mesh.C();
		d               = length(nelems);
		C.num_nodes     = prod([nelems...]*order + 1);
		C.num_elements  = prod([nelems...]);
		C.num_bdy_nodes = C.num_nodes - prod([nelems...]*order - 1);

		C.num_int_elements = prod([nelems...] - 1);
		C.num_bdy_elements = C.num_elements - C.num_int_elements;

		C.nnz = (order+2)^d*C.num_nodes;
		#       if (d == 2)
		#
		#       else
		#
		#       end
		return C
	end
	function ndgrid_fill(a, v, s, snext)
		for j = 1:length(a)
			a[j] = v[div(rem(j-1, snext), s)+1]
		end
	end
	function ndgrid{T}(vs::AbstractVector{T}...)
		n = length(vs)
		sz = map(length, vs)
		out = ntuple(n, i->Array(T, sz))
		s = 1
		for i=1:n
			a = out[i]::Array
			v = vs[i]
			snext = s*size(a,i)
			ndgrid_fill(a, v, s, snext)
			s = snext
		end
		out
	end
	function ismember(main_array, sub_array)
		out=int8(zeros(length(main_array)))
		match_index = findin(int64(main_array),int64(sub_array))
		out[match_index]=1
		out
	end
end
