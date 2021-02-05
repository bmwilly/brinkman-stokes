reload("stokes_flow/zero_bc.jl")

###NONZEROBC imposes Dirichlet boundary condition
function nonzerobc(A, f, xy, bound)
    nvtx = length(f); nbd = length(bound)
    null_col = spzeros(nvtx, nbd); null_row = spzeros(nbd, nvtx)
    Ax = A[1:nvtx, 1:nvtx];
    fx = f[1:nvtx]

  # set boundary condition
    xbd = xy[bound, 1]; ybd = xy[bound, 1]
    bc = zero_bc(xbd, ybd)
    fx -= Ax[:, bound] * bc
    dA = zeros(nvtx, 1); dA[bound] = ones(nbd, 2)
    Ax[:, bound] = null_col; Ax[bound, :] = null_row
    Ax += spdiagm(0 => vec(dA)); fx[bound] = bc
    return Ax, fx
end
