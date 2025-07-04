##SHAPE evaluates bilinear shape functions
#input
#    s           x coordinate
#    t           y coordinate
#output
#    phi         shape function
#    dphids      x derivative of phi
#    dphidt      y derivative of phi
function shape(s, t)

    phi = zeros(4)
    dphids = zeros(4)
    dphidt = zeros(4)

    phi[1] = 0.25 * (s - 1) * (t - 1)
    phi[2] = -0.25 * (s + 1) * (t - 1)
    phi[3] = 0.25 * (s + 1) * (t + 1)
    phi[4] = -0.25 * (s - 1) * (t + 1)

    dphids[1] = 0.25 * (t - 1)
    dphids[2] = -0.25 * (t - 1)
    dphids[3] = 0.25 * (t + 1)
    dphids[4] = -0.25 * (t + 1)

    dphidt[1] = 0.25 * (s - 1)
    dphidt[2] = -0.25 * (s + 1)
    dphidt[3] = 0.25 * (s + 1)
    dphidt[4] = -0.25 * (s - 1)

    (phi, dphids, dphidt)
end
