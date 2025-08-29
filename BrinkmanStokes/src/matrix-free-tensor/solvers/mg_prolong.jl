using SparseArrays

###MG_PROLONG GMG prolongation operator for square domain
# input
# 	nelx				number of elements in x-direction
# 	nely				number of elements in y-direction
# 	x 					x coordinate vector for coarse grid
# 	y 					y coordinate vector for coarse grid
# output
# 	P 					prolongation operator
function mg_prolong(nelx, nely, x, y)

    ecx = nelx / 2
    ecy = nely / 2

    tx = x[2:(nelx + 1)] - x[1:nelx]
    bx = zeros(int(nelx), 1)
    bx[1:2:(nelx - 1)] = x[3:2:(nelx + 1)] - x[1:2:(nelx - 1)]
    bx[2:2:nelx] = bx[1:2:(nelx - 1)]
    dx = tx ./ bx

    ty = y[2:(nely + 1)] - y[1:nely]
    by = zeros(int(nely), 1)
    by[1:2:(nely - 1)] = y[3:2:(nely + 1)] - y[1:2:(nely - 1)]
    by[2:2:nely] = by[1:2:(nely - 1)]
    dy = ty ./ by

    p = zeros(int(2ecx + 1), int(ecx + 1))
    # p = zeros(int(2ecx + 1), int(ecx + 1), int(3ecx + 1))
    # p = sparse(zeros(2ecx + 1), zeros(ecx + 1), 3ecx + 1)
    for j in 2:ecx
        p[2j - 2, j] = dx[2j - 3]
        p[2j - 1, j] = 1.0
        p[2j, j] = dx[2j]
    end
    j = 1
    p[2j - 1, j] = 1.0
    p[2j, 1] = dx[2j]
    j = ecx + 1
    p[2j - 2, j] = dx[2j - 3]
    p[2j - 1, j] = 1.0

    P = zeros(int((2ecy + 1) * (2ecx + 1)), int((ecy + 1) * (ecx + 1)))
    # P = zeros(int((2ecy + 1) * (2ecx + 1)), int((ecy + 1) * (ecx + 1)), int((3ecy + 1) * (3ecx + 1)))
    for j in 2:ecy
        cols = (j - 1) * (ecx + 1) + [1:(ecx + 1)]
        P[(2j - 3) * (2ecx + 1) + [1:(2ecx + 1)], cols] = dy[2j - 3] * p
        P[(2j - 2) * (2ecx + 1) + [1:(2ecx + 1)], cols] = p
        P[(2j - 1) * (2ecx + 1) + [1:(2ecx + 1)], cols] = dy[2j] * p
    end
    j = 1
    cols = (j - 1) * (ecx + 1) + [1:(ecx + 1)]
    P[(2j - 2) * (2ecx + 1) + [1:(2ecx + 1)], cols] = p
    P[(2j - 1) * (2ecx + 1) + [1:(2ecx + 1)], cols] = dy[2j] * p
    j = ecy + 1
    P[(2j - 3) * (2ecx + 1) + [1:(2ecx + 1)], cols] = dy[2j - 3] * p
    P[(2j - 2) * (2ecx + 1) + [1:(2ecx + 1)], cols] = p

    return P

end
