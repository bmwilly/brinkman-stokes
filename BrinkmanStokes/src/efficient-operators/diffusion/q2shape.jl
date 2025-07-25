###Q2SHAPE evaluates shape functions for a Q2 element
# input
#   r         x coordinate
#   s         y coordinate
# output
#   t         basis functions
#   dtdr      r basis derivatives
#   dtds      s basis derivatives
#
#  Reference Element Q9:
#
#    |
#    1  4--7--3
#    |  |     |
#    |  |     |
#    S  8  9  6
#    |  |     |
#    |  |     |
#    0  1--5--2
#    |
#    +--0--R--1-->
function q2shape(r, s)

    t = zeros(9, 1)
    dtdr = zeros(9, 1)
    dtds = zeros(9, 1)

    t[1] = 4.0 * (r - 1.0) * (r - 0.5) * (s - 1.0) * (s - 0.5)
    t[2] = 4.0 * r * (r - 0.5) * (s - 1.0) * (s - 0.5)
    t[3] = 4.0 * r * (r - 0.5) * s * (s - 0.5)
    t[4] = 4.0 * (r - 1.0) * (r - 0.5) * s * (s - 0.5)
    t[5] = -8.0 * r * (r - 1.0) * (s - 1.0) * (s - 0.5)
    t[6] = -8.0 * r * (r - 0.5) * s * (s - 1.0)
    t[7] = -8.0 * r * (r - 1.0) * s * (s - 0.5)
    t[8] = -8.0 * (r - 1.0) * (r - 0.5) * s * (s - 1.0)
    t[9] = 16.0 * r * (r - 1.0) * s * (s - 1.0)

    dtdr[1] = 4.0 * (2.0 * r - 1.5) * (s - 1.0) * (s - 0.5)
    dtdr[2] = 4.0 * (2.0 * r - 0.5) * (s - 1.0) * (s - 0.5)
    dtdr[3] = 4.0 * (2.0 * r - 0.5) * s * (s - 0.5)
    dtdr[4] = 4.0 * (2.0 * r - 1.5) * s * (s - 0.5)
    dtdr[5] = -8.0 * (2.0 * r - 1.0) * (s - 1.0) * (s - 0.5)
    dtdr[6] = -8.0 * (2.0 * r - 0.5) * s * (s - 1.0)
    dtdr[7] = -8.0 * (2.0 * r - 1.0) * s * (s - 0.5)
    dtdr[8] = -8.0 * (2.0 * r - 1.5) * s * (s - 1.0)
    dtdr[9] = 16.0 * (2.0 * r - 1.0) * s * (s - 1.0)

    dtds[1] = 4.0 * (r - 1.0) * (r - 0.5) * (2.0 * s - 1.5)
    dtds[2] = 4.0 * r * (r - 0.5) * (2.0 * s - 1.5)
    dtds[3] = 4.0 * r * (r - 0.5) * (2.0 * s - 0.5)
    dtds[4] = 4.0 * (r - 1.0) * (r - 0.5) * (2.0 * s - 0.5)
    dtds[5] = -8.0 * r * (r - 1.0) * (2.0 * s - 1.5)
    dtds[6] = -8.0 * r * (r - 0.5) * (2.0 * s - 1.0)
    dtds[7] = -8.0 * r * (r - 1.0) * (2.0 * s - 0.5)
    dtds[8] = -8.0 * (r - 1.0) * (r - 0.5) * (2.0 * s - 1.0)
    dtds[9] = 16.0 * r * (r - 1.0) * (2.0 * s - 1.0)

    t, dtdr, dtds

end
