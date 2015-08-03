using LinearOperators
reload("helpers/meshgrid.jl")
reload("helpers/input.jl")
reload("diffusion/deriv.jl")
reload("diffusion/qderiv.jl")
reload("diffusion/lderiv.jl")
reload("grids/gauleg.jl")
reload("grids/tprod.jl")
reload("grids/q2grid.jl")
reload("grids/q2grid_element_num.jl")
reload("grids/q2grid_node_num.jl")
reload("diffusion/q2deriv.jl")

# generate Q2 grid for square cavity

########################
# grid = cavity_domain()
########################
## define geometry
println("Grid generation for cavity domain.")
msize = 2
msize = int(input("Mesh size: "))
n = 2^msize
np = int(n/2)
nel = q2grid_element_num(np, np)

# y-direction
yy = [1/np:1/np:1]
ypos = [0, yy]
yneg = -yy[length(yy):-1:1]
y = [yneg, ypos]'
x = y

# compute biquadratic element coordinates
nvtx = q2grid_node_num(np, np)
(X,Y) = meshgrid(x,y)
xx = reshape(X', nvtx, 1)
yy = reshape(Y', nvtx, 1)
xy = [xx[:] yy[:]]

mv = q2grid(np, np)


##############################
# grid = q2p1grid(cavity_grid)
##############################
grid = {
  "x" => x,
  "y" => y,
  "xy" => xy,
  "mv" => mv
}

# stokes q2-p1 matrix generator

########################################
# stokes_mats = stokes_q2p1(stokes_grid)
########################################
xy = grid["xy"];
mv = grid["mv"]

nngpt = 9
x = xy[:, 1]
y = xy[:, 2]
nvtx = length(x)
nu = 2nvtx
nel = length(mv[:, 1])

println("setting up Q2 matrices... ")

# initialize global matrices
A = spzeros(nu, nu)

## Gauss point integration rules
# 3x3 Gauss points
quadrule_1D = gauleg(-1, 1, 3)
quadrule_2D = tprod(quadrule_1D)
s = quadrule_2D["x"][:,1]
t = quadrule_2D["x"][:,2]
wt = quadrule_2D["w"]

# inner loop over elements
xlv = zeros(nel, 4)
ylv = zeros(nel, 4)
for ivtx = 1:4
    xlv[:, ivtx] = x[mv[:, ivtx]]
    ylv[:, ivtx] = y[mv[:, ivtx]]
end

# initialize derivative matrices
ae = zeros(nel, 9, 9)

# loop over Gauss points
for igpt = 1:nngpt
    sigpt = s[igpt]
    tigpt = t[igpt]
    wght = wt[igpt]

    # evaluate derivatives, etc.
    (jac, invjac, phi, dphidx, dphidy) = deriv(sigpt, tigpt, xlv, ylv)
    (psi, dpsidx, dpsidy) = q2deriv(sigpt, tigpt, xlv, ylv)
    for j = 1:9
        for i = 1:9
            ae[:, i, j] += wght * dpsidx[:, i] .* dpsidx[:, j] .* invjac[:]
            ae[:, i, j] += wght * dpsidy[:, i] .* dpsidy[:, j] .* invjac[:]
        end
    end
end # end of Gauss point loop

## element assembly into global matrices
# component velocity matrices
for krow = 1:9
  nrow = mv[:, krow]
  for kcol = 1:9
    ncol = mv[:, kcol]
    A += sparse(nrow, ncol, ae[:, krow, kcol], nu, nu)
    A += sparse(nrow + nvtx, ncol + nvtx, ae[:, krow, kcol], nu, nu)
  end
end

a_density = nnz(A) / (size(A, 1) * size(A, 2))
println("A density: $(a_density)")

@time for cnt = 1:100; u = rand(size(A,1), 1); w = A*u; end

# nza = float(bool(A))
# nzb = float(bool(u))
# f = nza * nzb
# f = 2f - float(bool(f))
# f = sum(sum(f))
#
# println("# flops: $(f)")
