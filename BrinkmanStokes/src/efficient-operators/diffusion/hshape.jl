function hshape(r, s, p)

    nn = (p + 1) * (p + 1)
    t = rand(nn, 1)
    dtdr = rand(nn, 1)
    dtds = rand(nn, 1)

    t, dtdr, dtds
end
