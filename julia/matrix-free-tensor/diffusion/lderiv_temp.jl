include("lshape.jl")

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

# evaluate bilinear shape functions
(phie, dphids, dphidt) = shape(s, t)

# evaluate linear shape functions
(chie, dchids, dchidt) = lshape(s, t)

# local derivatives
dxds = copy(zerov)
dxdt = copy(zerov)
dyds = copy(zerov)
dydt = copy(zerov)
invjac = copy(zerov)

for ivtx = 1:4
    dxds[:] += xl[:,ivtx] .* onev * dphids[ivtx]
    dxdt[:] += xl[:,ivtx] .* onev * dphidt[ivtx]
    dyds[:] += yl[:,ivtx] .* onev * dphids[ivtx]
    dydt[:] += yl[:,ivtx] .* onev * dphidt[ivtx]
end

chi = zeros(nel, 3)
dchidx = zeros(nel, 3)
dchidy = zeros(nel, 3)

for ivtx = 1:3
    chi[:, ivtx] = chie[ivtx] * onev
    dchidx[:, ivtx] = dchids[ivtx] .* dydt[:] - dchidt[ivtx] .* dyds[:]
    dchidy[:, ivtx] = -dchids[ivtx] .* dxdt[:] + dchidt[ivtx] .* dxds[:]
end

(chi, dchidx, dchidy)
