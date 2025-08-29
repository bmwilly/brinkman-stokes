###MG_ZEROBC imposes zero boundary condition
# input
# 	A 				matrix with no boundary conditions imposed
# 	xy
# 	bound 			index set of Dirichlet boundary nodes
# output
# 	Ax 				modified matrix with boundary conditions imposed
function mg_zerobc(A, xy, bound)

    nvtx = length(A[1, :]); nbd = length(bound)
    nullcol = spzeros(nvtx, nbd); nullrow = spzeros(nbd, nvtx)
    Ax = A[1:nvtx, 1:nvtx]

    # impose boundary condition
    dA = zeros(nvtx, 1)
    dA[bound] = ones(nbd, 1)
    Ax[:, bound] = nullcol
    Ax[bound, :] = nullrow
    return Ax += spdiagm(0 => vec(dA))

end
