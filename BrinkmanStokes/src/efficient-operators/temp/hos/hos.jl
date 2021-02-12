using LinearOperators
reload("helpers/meshgrid.jl")
reload("helpers/input.jl")
reload("grids/gauleg.jl")
reload("grids/tprod.jl")
reload("grids/qgrid.jl")
reload("grids/grid_element_num.jl")
reload("grids/grid_node_num.jl")
reload("diffusion/hderiv.jl")
reload("diffusion/deriv.jl")
reload("diffusion/qderiv.jl")
reload("diffusion/lderiv.jl")
reload("diffusion/gauss_source.jl")
reload("diffusion/nonzerobc.jl")
reload("diffusion/hosbc.jl")

########################
# grid = cavity_domain()
########################
## define geometry
println("Grid generation for cavity domain.")
# msize = 2
p = int(input("Polynomial order: "))
msize = int(input("Mesh size: "))
n = 2^msize
np = int(n / 2)
nel = grid_element_num(np, np)

# y-direction
# yy = [(p-1)/(p*np):(p-1)/(p*np):1]
yy = [2 / (p * np):2 / (p * np):1]
ypos = [0, yy]
yneg = -yy[length(yy):-1:1]
y = [yneg, ypos]'
x = y

# compute bicubic element coordinates
nvtx = grid_node_num(np, np, p)
(X, Y) = meshgrid(x, y)
xx = reshape(X', nvtx, 1)
yy = reshape(Y', nvtx, 1)
xy = [xx[:] yy[:]]

mv = qgrid(np, np, p)

# compute boundary vertices and edges
# four boundary edges
k1 = findall(xy[:,2] .== -1)
k2 = findall((xy[:,1] .== 1) .& (xy[:,2] .< 1) .& (xy[:,2] .> -1))
k3 = findall(xy[:,2] .== 1)
k4 = findall((xy[:,1] .== -1) .& (xy[:,2] .< 1) .& (xy[:,2] .> -1))
bound = sort([k1; k2; k3; k4])

##############################
# grid = q2p1grid(cavity_grid)
##############################

# stokes q2-p1 matrix generator

########################################
# stokes_mats = stokes_q2p1(stokes_grid)
########################################
nngpt = (p + 1) * (p + 1)
x = xy[:, 1]
y = xy[:, 2]
nvtx = length(x)
nu = 2nvtx
nel = length(mv[:, 1])

# initialize global matrices
A = spzeros(nu, nu)
f = zeros(nu, 1)

## Gauss point integration rules
# (p+1)x(p+1) Gauss points
quadrule_1D = gauleg(-1, 1, p + 1)
quadrule_2D = tprod(quadrule_1D)
s = quadrule_2D["x"][:,1]
t = quadrule_2D["x"][:,2]
wt = quadrule_2D["w"]

# inner loop over elements
xlv = zeros(nel, p * p)
ylv = zeros(nel, p * p)
for ivtx = 1:(p * p)
    xlv[:, ivtx] = x[mv[:, ivtx]]
    ylv[:, ivtx] = y[mv[:, ivtx]]
end

# initialize derivative matrices
ae = zeros(nel, nngpt, nngpt)
# fe = zeros(nel, 9)

# loop over Gauss points
for igpt = 1:nngpt
    sigpt = s[igpt]
    tigpt = t[igpt]
    wght = wt[igpt]

    # evaluate derivatives, etc.
    (jac, invjac, phi, dphidx, dphidy) = deriv(sigpt, tigpt, xlv, ylv)
    (psi, dpsidx, dpsidy) = qderiv(sigpt, tigpt, xlv, ylv)
    (chi, dchidx, dchidy) = lderiv(sigpt, tigpt, xlv, ylv)
    # rhs = gauss_source(sigpt, tigpt, xlv, ylv)
    # (jac, invjac, chi, dchidx, dchidy) = hderiv(sigpt, tigpt, xlv, ylv)
    for j = 1:nngpt
        for i = 1:nngpt
            ae[:, i, j] += wght * dpsidx[:, i] .* dpsidx[:, j] .* invjac[:]
            ae[:, i, j] += wght * dpsidy[:, i] .* dpsidy[:, j] .* invjac[:]
        end
        # fe[:, j] += wght * rhs[:] .* psi[:, j] .* jac[:]
    end
end # end of Gauss point loop

## element assembly into global matrices
# component velocity matrices
@time (for krow = 1:nngpt
    nrow = mv[:, krow]
    for kcol = 1:nngpt
        ncol = mv[:, kcol]
        A += sparse(nrow, ncol, ae[:, krow, kcol], nu, nu)
        A += sparse(nrow .+ nvtx, ncol .+ nvtx, ae[:, krow, kcol], nu, nu)
    end
  # f[nrow, 1] += fe[:, krow]
end)

# Agal,fgal = nonzerobc(A, f, xy, bound)
Ast, fst = hosbc(A, f, xy, bound)

a_density = nnz(Ast) / (size(Ast, 1) * size(Ast, 2))
println("A density: $(a_density)")

@time for cnt = 1:100; u = rand(size(Ast, 1), 1); w = Ast * u; end

# Ax = A[1:nvtx, 1:nvtx]; Ay = A[nvtx + 1:nu, nvtx + 1:nu]
# u = rand(size(Ast, 1), 1)
# ux = u[1:nvtx]; uy = u[nvtx+1:nu]
# wx = Ax*ux; wy = A

# nza = float(bool(A))
# nzb = float(bool(u))
# f = nza * nzb
# f = 2f - float(bool(f))
# f = sum(sum(f))
#
# println("# flops: $(f)")
