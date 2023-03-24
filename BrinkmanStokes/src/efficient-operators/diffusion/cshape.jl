###CSHAPE evaluates biquadratic shape functions
#input
#    s           x coordinate
#    t           y coordinate
#output
#    chi         shape function
#    dchids      x derivative of chi
#    dchidt      y derivative of chi
function cshape(s, t)

    ellx = zeros(3, 1)
    elly = zeros(3, 1)
    dellx = zeros(3, 1)
    delly = zeros(3, 1)

    psi = zeros(9, 1)
    dpsids = zeros(9, 1)
    dpsidt = zeros(9, 1)

    # one-dimensional shape functions
    ellx[1] = 0.5s * (s - 1)
    elly[1] = 0.5t * (t - 1)
    ellx[2] = 1 - (s * s)
    elly[2] = 1 - (t * t)
    ellx[3] = 0.5s * (s + 1)
    elly[3] = 0.5t * (t + 1)

    dellx[1] = s - 0.5
    delly[1] = t - 0.5
    dellx[2] = -2s
    delly[2] = -2t
    dellx[3] = s + 0.5
    delly[3] = t + 0.5

    # two-dimensional shape functions
    psi[1] = ellx[1] * elly[1]
    psi[2] = ellx[3] * elly[1]
    psi[3] = ellx[3] * elly[3]
    psi[4] = ellx[1] * elly[3]
    psi[5] = ellx[2] * elly[1]
    psi[6] = ellx[3] * elly[2]
    psi[7] = ellx[2] * elly[3]
    psi[8] = ellx[1] * elly[2]
    psi[9] = ellx[2] * elly[2]

    dpsids[1] = dellx[1] * elly[1]
    dpsids[2] = dellx[3] * elly[1]
    dpsids[3] = dellx[3] * elly[3]
    dpsids[4] = dellx[1] * elly[3]
    dpsids[5] = dellx[2] * elly[1]
    dpsids[6] = dellx[3] * elly[2]
    dpsids[7] = dellx[2] * elly[3]
    dpsids[8] = dellx[1] * elly[2]
    dpsids[9] = dellx[2] * elly[2]

    dpsidt[1] = ellx[1] * delly[1]
    dpsidt[2] = ellx[3] * delly[1]
    dpsidt[3] = ellx[3] * delly[3]
    dpsidt[4] = ellx[1] * delly[3]
    dpsidt[5] = ellx[2] * delly[1]
    dpsidt[6] = ellx[3] * delly[2]
    dpsidt[7] = ellx[2] * delly[3]
    dpsidt[8] = ellx[1] * delly[2]
    dpsidt[9] = ellx[2] * delly[2]

    # three-dimensional shape function

    (psi, dpsids, dpsidt)
end
