reload("../julia-homg/Basis.jl")
reload("../julia-homg/Refel.jl")

function hshape(r, s, p)

    nn = (p + 1) * (p + 1)
    t = zeros(nn, 1)
    dtdr = zeros(nn, 1)
    dtds = zeros(nn, 1)

    t = rand(nn, 1)
    dtdr = rand(nn, 1)
    dtds = rand(nn, 1)


    return t, dtdr, dtds
end
