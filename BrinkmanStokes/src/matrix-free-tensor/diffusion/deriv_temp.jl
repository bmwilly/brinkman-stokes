include("shape.jl")

igpt = 1
sigpt = s[igpt]
tigpt = t[igpt]
wght = wt[igpt]

xl = xlv
yl = ylv
s = sigpt
t = tigpt

nel = length(xl[:, 1])
zerov = zeros(nel, 1)
onev = ones(nel, 1)

# evaluate shape functions
(phie, dphids, dphidt) = shape(s, t)

dxds = copy(zerov)
dxdt = copy(zerov)
dyds = copy(zerov)
dydt = copy(zerov)
jac = copy(zerov)
invjac = copy(zerov)

for ivtx = 1:4
    dxds[:] += xl[:, ivtx] .* onev * dphids[ivtx]
    dxdt[:] += xl[:, ivtx] .* onev * dphidt[ivtx]
    dyds[:] += yl[:, ivtx] .* onev * dphids[ivtx]
    dydt[:] += yl[:, ivtx] .* onev * dphidt[ivtx]
end

jac[:] = dxds[:].*dydt[:] - dxdt[:].*dyds[:]

# check element jacobian
if any(jac .< 1e-9)
    println("Bad element warning...")
    if any(jac .<= 0.0)
        error("singular Jacobian ... aborted ...")
    end
end

invjac[:] = onev ./ jac[:]

phi = zeros(nel, 4)
dphidx = zeros(nel, 4)
dphidy = zeros(nel, 4)

for ivtx = 1:4
    phi[:, ivtx] = phie[ivtx] * onev
    dphidx[:, ivtx] = dphids[ivtx] .* dydt[:] - dphidt[ivtx] .* dyds[:]
    dphidy[:, ivtx] = -dphids[ivtx] .* dxdt[:] + dphidt[ivtx] .* dxds[:]
end

(jac, invjac, phi, dphidx, dphidy)
