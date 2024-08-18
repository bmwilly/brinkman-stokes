reload("diffusion/unit_rhs.jl")

###GAUSS_SOURCE evaluates source term at Gauss point
function gauss_source(s, t, xl, yl)
    nel = length(xl[:, 1])
    xx = zeros(nel, 1)
    yy = zeros(nel, 1)
    (phie, dphids, dphidt) = shape(s, t)
    for ivtx = 1:4
        xx += phie[ivtx] * xl[:, ivtx]
        yy += phie[ivtx] * yl[:, ivtx]
    end
    ff = unit_rhs(xx, yy, nel)
end
