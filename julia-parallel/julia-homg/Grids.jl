
module Grids
include("Refel.jl")
include("Basis.jl")
include("Tensor.jl")
include("BaseCustom.jl")
include("Hexmesh.jl")

type Grid
	level
    is_finest
    eig_max
    eig_min

    k_evec

    jacobi_omega
    jacobi_invdiag
    jacobi_inv_block_diag

    ssor_M
    ssor_N
    sor_omega

    smoother
    linear_smoother
    K
    L
    K_lin
    Boundary

    M

    R
    P

    Mesh

    Coarse  # handle to coarse grid

	function Grid(mesh, order, coarse=[])
		grid = new()
		if typeof(coarse) != Grid
			grid.level = 0;
			grid.Coarse = [];
		else
			grid.level = coarse.level + 1;
			grid.Coarse = coarse;
		end

		grid.Mesh = mesh;

		Mesh.set_order(mesh, order);
		grid.sor_omega = 1;
		if typeof(coarse) == Grid
			grid.P = Mesh.assemble_interpolation(grid.Coarse.Mesh, order);
			grid.R = grid.P';
		end

		# for Dirichlet boundary conditions
		grid.Boundary = Mesh.get_boundary_node_indices(mesh, order);

		## defaults ...
		grid.smoother = "sor";
		grid.jacobi_omega = 2/3;

		grid.jacobi_invdiag = [];
		grid.R = [];
		grid.eig_max = [];
		grid.ssor_M = [];
		grid.ssor_N = [];

		return grid
	end
end
function assemble_poisson(grid, mu)
	# fine grid material props ...
	Mesh.set_coeff (grid.Mesh, mu) ;
	# assemble for this level ...
	(grid.K, grid.M, grid.jacobi_inv_block_diag) = Mesh.assemble_poisson(grid.Mesh, grid.Mesh.order);
	# syms x y z
	if ( grid.Mesh.dim == 2 )
		fx = (x,y)->(-8*pi^2*(sin(2*pi*x) .* sin(2*pi*y)))
	else
		fx = (x,y)->(-12*pi^2*(sin(2*pi*x) .* sin(2*pi*y) * sin(2*pi*z) ))
	end

	grid.L = Mesh.assemble_rhs(grid.Mesh, fx, grid.Mesh.order);
	grid.L[grid.Boundary] = 0;

	# propagate to lower grids
	if typeof(grid.Coarse) == Grid
		if typeof(mu) != Function && isinteger(mu)
			harmonic = 0;
			if (grid.Mesh.dim == 2)
				mu2 = reshape(mu, grid.Mesh.nelems(1), grid.Mesh.nelems(2));
				if (harmonic != 0) #if (harmonic)
					mu_coarse = 4 ./ ( 1./mu2[1:2:end, 1:2:end] + 1./mu2[2:2:end, 1:2:end] + 1./mu2[1:2:end, 2:2:end] + 1./mu2[2:2:end, 2:2:end] );
				else
					mu_coarse = 0.25*(mu2[1:2:end, 1:2:end] + mu2[2:2:end, 1:2:end] + mu2[1:2:end, 2:2:end] + mu2[2:2:end, 2:2:end]);
				end
			else
				mu3 = reshape(mu, grid.Mesh.nelems(1), grid.Mesh.nelems(2), grid.Mesh.nelems(3));
				if (harmonic != 0) #if (harmonic)
					mu_coarse = 8 ./ ( 1./mu3[1:2:end, 1:2:end, 1:2:end] + 1./mu3[2:2:end, 1:2:end, 1:2:end] +
						1./mu3[1:2:end, 2:2:end, 1:2:end] + 1./mu3[2:2:end, 2:2:end, 1:2:end] +
						1./mu3[1:2:end, 1:2:end, 2:2:end] + 1./mu3[2:2:end, 1:2:end, 2:2:end] +
						1./mu3[1:2:end, 2:2:end, 2:2:end] + 1./mu3[2:2:end, 2:2:end, 2:2:end] );
				else
					mu_coarse = 0.125*(mu3[1:2:end, 1:2:end, 1:2:end] + mu3[2:2:end, 1:2:end, 1:2:end] + mu3[1:2:end, 2:2:end, 1:2:end] + mu3[2:2:end, 2:2:end, 1:2:end] +
						mu3[1:2:end, 1:2:end, 2:2:end] + mu3[2:2:end, 1:2:end, 2:2:end] + mu3[1:2:end, 2:2:end, 2:2:end] + mu3[2:2:end, 2:2:end, 2:2:end] );
				end
			end
			assemble_poisson(grid.Coarse, mu_coarse[:]) ;
		else
			assemble_poisson(grid.Coarse, mu) ;
		end
	end
end

function use_linearized_smoothers(grid)
	grid.linear_smoother = true;
	(grid.K_lin, M)  =  Mesh.assemble_poisson_linearized(grid.Mesh, grid.Mesh.order);
	if (typeof(grid.Coarse) == Grid )
		use_linearized_smoothers(grid.Coarse);
	end
end

function set_stiffness(grid, K)
	grid.K = K;
end

# compute the residual
function residual(grid, rhs, u)

	if typeof(rhs) == Nothing && typeof(u) == Nothing
		rhs = grid.L;
	end
	if typeof(u) == Nothing
		u = zeros(size(rhs));
	end

	r = grid.K*u - rhs;
	return r

end

function solve_pcg(grid, num_vcyc, smoother, v1, v2, rhs, u)
	set_smoother(grid, smoother);

	r = residual(grid, rhs, u);
	rho = zeros(size(u));

	rho = vcycle(grid, v1, v2, r, rho);
	p = rho;
	display(string("Initial residual is ", norm(r)));
	display("------------------------------------------");
	r0 = norm(r);
	for i=1:num_vcyc
		h = grid.K * p;
		rho_res = dot(rho[:], r[:]);
		alpha = rho_res / dot(p[:], h[:] );
		u = u + alpha*p;
		r = r - alpha*h;

		display(string(i, ": |res| = ", norm(r)));
		if (norm(r)/r0 < 1e-8)
			iter = i;
			rr = norm(r)/r0;
			return;
		end

		# precondition ..
		rho = zeros(size(u));
		rho = vcycle(grid, v1, v2, r, rho);

		beta = dot(rho[:], r[:]) / rho_res ;
		p = rho + beta*p;
	end
	display("------------------------------------------");
	iter = num_vcyc;
	rr = norm(r)/r0;
	return u, rr, iter
end

function solve(grid, num_vcyc, smoother, v1, v2, rhs, u)
	set_smoother(grid, smoother);
	r = residual(grid, rhs, u);
	display(string("Initial residual is ", norm(r)));
	display("------------------------------------------");
	r0 = norm(r);
	for i=1:num_vcyc
		u = vcycle(grid, v1, v2, rhs, u);
		r = residual(grid, rhs, u);
		display(string(i, ": |res| = ", norm(r)));
		if (norm(r)/r0 < 1e-8)
			iter = i;
			rr = norm(r)/r0;
			return;
		end
	end
	display("------------------------------------------");
	iter = num_vcyc;
	rr = norm(r)/r0;
	return u, rr, iter
end

function vcycle(grid, v1, v2, rhs, u)
	# function u = vcycle(grid, v1, v2, rhs, u)
	# solve system using initial guess u, given rhs
	# with v1 pre and v2 post-smoothing steps

	# handle for the coarsest level
	if is(grid.Coarse, nothing)
		u = grid.K \ rhs;
		return;
	end
	# 1. pre-smooth
	u = smooth(grid, v1, rhs, u );

	# 2. compute residual
	res = residual(grid, rhs, u);

	# 3. restrict
	if ~isempty(grid.R)
		res_coarse = grid.R * res;
		res_coarse[grid.Coarse.Boundary] = 0;

		# 4. recurse
		u_corr_coarse = vcycle(grid.Coarse, v1, v2, res_coarse, zeros(size(res_coarse)));

		# 5. prolong and correct
		u = u - grid.P * u_corr_coarse;
	end
	# 6. post-smooth
	u = smooth(grid, v2, rhs, u);
	return u;
end # v-cycle


# smoothers
function smooth(grid, v, rhs, u)
	if grid.smoother == "jacobi"
		u = smoother_jacobi(grid, v, rhs, u);
		return u;
	elseif grid.smoother == "blk_jac"
		u = smoother_block_jacobi(grid, v, rhs, u);
		return u;
	elseif grid.smoother == "chebyshev"
		u = smoother_chebyshev(grid, v, rhs, u);
		return u;
	elseif grid.smoother == "ssor"
		u = smoother_sym_sor(grid, v, rhs, u);
		return u;
	else
		display("ERROR: Unrecognized smoother type");
		return;
	end
end

function set_coeff(grid, mu)
	Mesh.set_coeff (grid.Mesh, mu) ;
	if (typeof(grid.Coarse) == Grid )
		Mesh.set_coeff (grid.Coarse.Mesh, mu) ;
	end
end

function set_smoother(grid, sm)
	grid.smoother = sm;
	if (typeof(grid.Coarse) == Grid )
		set_smoother(grid.Coarse, sm);
	end
end

function smoother_jacobi (grid, v, rhs, u)
	# standard jacobi smoother
	if isempty(grid.jacobi_invdiag)
		D = diag(grid.K);
		grid.jacobi_invdiag = 1./D;
	end

	for i=1:v
		res  = grid.jacobi_invdiag .* residual(grid, rhs, u);
		u = u - grid.jacobi_omega.*res;
	end
	return u
end # jacobi

function smoother_block_jacobi (grid, v, rhs, u)
	# block jacobi smoother
	if ( isempty(grid.jacobi_inv_block_diag) )
		error("inv block doagonal not assembled");
	end

	for i=1:v
		res  = grid.jacobi_inv_block_diag * residual(grid, rhs, u);
		u = u - grid.jacobi_omega.*res;
	end
	return u
end # blk jacobi

function smoother_sym_sor (grid, v, rhs, u)
	if ( isempty ( grid.ssor_M ) )
		w = grid.sor_omega;
		n = length(u);
		grid.ssor_M = spdiagm( (1/w)*diag(grid.K), 0, n, n) + tril(grid.K,-1);
		grid.ssor_N = spdiagm(((1-w)/w)*diag(grid.K), 0, n, n) - triu(grid.K,1);
	end

	for i=1:v
		r = residual(grid, rhs, u);
		u = u - grid.ssor_M \ r;
		u = grid.ssor_M' \ (grid.ssor_N'*u + rhs);
	end
	return u
end

function set_sor_omega(grid, w)
	grid.sor_omega = w;
	grid.ssor_M = [];
	grid.ssor_N = [];
	return grid
end

function smoother_chebyshev (grid, v, rhs, u)
	if ( isempty ( grid.eig_max ) )
		D = diag(grid.K);
		grid.jacobi_invdiag = 1./D;
		Kc = spdiagm(grid.jacobi_invdiag,[0],length(D), length(D)) * grid.K;
		opts.tol = 0.01;
		grid.eig_max = eigs(Kc, 1, which="LM", opts);
		# grid.eig_min = eigs(Kc, 1, 'sm');
	end

	# adjust the eigenvalues to hit the upper spectrum
	beta = grid.eig_max;
	alpha = 0.25*grid.eig_max;# (grid.eig_min + grid.eig_max)/2;

	delta = (beta - alpha)/2;
	theta = (beta + alpha)/2;
	s1 = theta/delta;
	rhok = 1./s1;

	d = zeros(size(u));

	# first loop
	res = -residual (grid, rhs, u );
	d = res/theta.* grid.jacobi_invdiag;
	u = u + d;

	for iter = 2:v
		rhokp1 = 1/ (2*s1 - rhok);
		d1 = rhokp1 * rhok;
		d2 = 2*rhokp1 / delta;
		rhok = rhokp1;
		res = -residual (grid, rhs, u );
		d = d1 * d + d2 * res.*grid.jacobi_invdiag;
		u = u + d;
	end
	return u
end # chebyshev

function get_eigenvectors(grid)
  # generate the correct matrix
  Kc = grid.K; #(eye(size(grid.K)) - grid.ZeroBoundary) + grid.ZeroBoundary * grid.K * grid.ZeroBoundary;
  (evec, eval) = svd(full(Kc)); #eig(full(Kc), full(grid.M));
  (eval,per)= sort(diag(eval)); #'ascend'
  evec = evec[:,per];
  return evec, eval
end


function get_u0(grid)
  u0 = rand(size(grid.L));
  u0[grid.Boundary] = 0;
  return u0;
end

#=
function spdiagm{T<:Integer}(B, d::Array{T, 1}, m::Integer, n::Integer)
	ndiags = length(d)
	ncoeffs = 0
	for vec in B
		ncoeffs += length(vec)
	end
	I = Array(T, ncoeffs)
	J = Array(T, ncoeffs)
	V = Array(Float64, ncoeffs)
	id = 0
	i = 0
	for vec in B
		id += 1
		diag = d[id]
		numel = length(vec)
		if diag < 0
			row = -diag
			col = 0
		elseif diag > 0
			row = 0
			col = diag
		else
			row = 0
			col = 0
		end
		range = 1+i:numel+i
		I[range] = row+1:row+numel
		J[range] = col+1:col+numel
		V[range] = vec[1:numel]
		i += numel
	end
	return sparse(I, J, V, m, n)
end
=#

end
