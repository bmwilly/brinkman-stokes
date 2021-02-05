using LinearOperators
reload("helpers/meshgrid.jl")
reload("helpers/input.jl")
reload("diffusion/deriv.jl")
reload("diffusion/qderiv.jl")
reload("diffusion/lderiv.jl")

# generate Q2 grid for square cavity

########################
# grid = cavity_domain()
########################
## define geometry
println("Grid generation for cavity domain.")
msize = 2
msize = int(input("Mesh size: "))
n = 2^msize
np = int(n / 2)
nel = int(np^2)

# y-direction
yy = [1 / np:1 / np:1]
ypos = [0, yy]
yneg = -yy[length(yy):-1:1]
y = [yneg, ypos]'
x = y

# compute biquadratic element coordinates
nvtx = (n + 1) * (n + 1)
(X, Y) = meshgrid(x, y)
xx = reshape(X', nvtx, 1)
yy = reshape(Y', nvtx, 1)
xy = [xx[:] yy[:]]

kx = 1
ky = 1
mel = 0
mv = zeros(Int64, nel, 9)
for j = 1:np
    for i = 1:np
        mref = (n + 1) * (ky - 1) + kx
        println(mref)
        mel += 1
        nvv = zeros(9)
        nvv[1] = mref
        nvv[2] = mref + 2
        nvv[3] = mref + 2n + 4
        nvv[4] = mref + 2n + 2
        nvv[5] = mref + 1
        nvv[6] = mref + n + 3
        nvv[7] = mref + 2n + 3
        nvv[8] = mref + n + 1
        nvv[9] = mref + n + 2
        mv[mel, 1:9] = nvv[1:9]
        kx += 2
    end
    ky += 2
    kx = 1
end

# compute boundary vertices and edges
# four boundary edges
k1 = findall(xy[:,2] .== -1)
e1 = Int[]
for k = 1:mel
    if any(mv[k,5] .== k1)
        push!(e1, k)
    end
end
ef1 = ones(size(e1))

k2 = findall((xy[:,1] .== 1) & (xy[:,2] .< 1) & (xy[:,2] .> -1))
e2 = Int[]
for k = 1:mel
    if any(mv[k,6] .== k2)
        push!(e2, k)
    end
end
ef2 = 2 * ones(size(e2))

k3 = findall(xy[:,2] .== 1)
e3 = Int[]
for k = 1:mel
    if any(mv[k,7] .== k3)
        push!(e3, k)
    end
end
ef3 = 3 * ones(size(e3))

k4 = findall((xy[:,1] .== -1) & (xy[:,2] .< 1) & (xy[:,2] .> -1))
e4 = Int[]
for k = 1:mel
    if any(mv[k,8] .== k4)
        push!(e4, k)
    end
end
ef4 = 4 * ones(size(e4))

bound = sort([k1; k2; k3; k4])
mbound = [e1' ef1'; e2' ef2'; e3' ef3'; e4' ef4']

grid = {
  "mv" => mv,
  "xy" => xy,
  "bound" => bound,
  "mbound" => mbound,
  "x" => x,
  "y" => y
}

##############################
# grid = q2p1grid(cavity_grid)
##############################
x = grid["x"]; y = grid["y"]; xy = grid["xy"];
mv = grid["mv"]; bound = grid["bound"]

## centroid coordinate vector
xx = xy[:, 1]
yy = xy[:, 2]
nvtx = length(xx)
nel = length(mv[:, 1])

## recompute mid-side points in the case of stretched grids
# y-direction
yv = yy
ny = length(y)

for k = 2:2:ny
    yold = y[k]
    ynew = 0.5 * (y[k + 1] + y[k - 1])
    l = findall(yy == yold)
    yv[l] .= ynew
    y[k] = ynew
end

# x-direction
xv = xx
nx = length(x)

for k = 2:2:nx
    xold = x[k]
    xnew = 0.5 * (x[k + 1] + x[k - 1])
    l = findall(xx == xold)
    xv[l] .= xnew
    x[k] = xnew
end

xy = [xv yv]

# centroid coordinates
xc = zeros(nel, 1)
yc = zeros(nel, 1)
for ielem = 1:nel
    xc[ielem] = mean(xx[mv[ielem, 1:4]])
    yc[ielem] = mean(yy[mv[ielem, 1:4]])
end

xyp = [xc yc]

# compute edge to edge connection array ee
np = nel

grid = {
  "x" => x,
  "y" => y,
  "xy" => xy,
  "xyp" => xyp,
  "mv" => mv
}

# stokes q2-p1 matrix generator

########################################
# stokes_mats = stokes_q2p1(stokes_grid)
########################################
xy = grid["xy"]; xyp = grid["xyp"]; mv = grid["mv"]

nngpt = 9
x = xy[:, 1]
y = xy[:, 2]
xp = xyp[:, 1]
yp = xyp[:, 2]
nvtx = length(x)
nu = 2nvtx
np = 3length(xp)
nel = length(mv[:, 1])
mp = [[1:3:3nel] [2:3:3nel] [3:3:3nel]]

println("setting up Q2-P1 matrices... ")

# initialize global matrices
A = spzeros(nu, nu)

## Gauss point integration rules
# 3x3 Gauss points
s = zeros(nngpt, 1)
t = zeros(nngpt, 1)
wt = zeros(nngpt, 1)
gpt = sqrt(0.6);
s[1] = -gpt; t[1] = -gpt; wt[1] = 25 / 81;
s[2] =  gpt; t[2] = -gpt; wt[2] = 25 / 81;
s[3] =  gpt; t[3] =  gpt; wt[3] = 25 / 81;
s[4] = -gpt; t[4] =  gpt; wt[4] = 25 / 81;
s[5] =  0.0; t[5] = -gpt; wt[5] = 40 / 81;
s[6] =  gpt; t[6] =  0.0; wt[6] = 40 / 81;
s[7] =  0.0; t[7] =  gpt; wt[7] = 40 / 81;
s[8] = -gpt; t[8] =  0.0; wt[8] = 40 / 81;
s[9] =  0.0; t[9] =  0.0; wt[9] = 64 / 81;

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
    (psi, dpsidx, dpsidy) = qderiv(sigpt, tigpt, xlv, ylv)
    (chi, dchidx, dchidy) = lderiv(sigpt, tigpt, xlv, ylv)
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
        A += sparse(nrow .+ nvtx, ncol .+ nvtx, ae[:, krow, kcol], nu, nu)
    end
end

mats = {
  "A" => A
}

mats = merge(mats, grid)

u = linspace(1, nu, nu)
A * u
