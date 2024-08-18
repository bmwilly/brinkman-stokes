###LSHAPE evaluates linear shape functions
#input
#    s           x coordinate
#    t           y coordinate
#output
#    chi         shape function
#    dchids      x derivative of chi
#    dchidt      y derivative of chi
function lshape(s, t)

    chi = zeros(3, 1)
    dchids = zeros(3, 1)
    dchidt = zeros(3, 1)

    chi[1] = 1.0
    chi[2] = s
    chi[3] = t

    dchids[1] = 0.0
    dchids[2] = 1.0
    dchids[3] = 0.0

    dchidt[1] = 0.0
    dchidt[2] = 0.0
    dchidt[3] = 1.0

    (chi, dchids, dchidt)
end
