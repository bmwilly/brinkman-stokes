reload("stokes_flow/zero_bc.jl")
reload("stokes_flow/obstacle_bc.jl")
reload("stokes_flow/poiseuille_bc.jl")

###STREAMBC imposes Dirichlet BC on the stream function
# input
#   A       stiffness matrix
#   f       RHS vector
#   xy      vertex coordinate vector
#   bound   boundary vertex vector
# output
#   Agal    stiffness matrix
#   fgal    RHS vector
function streambc(a, f, xy, bound, domain)

  nvtx = length(f); nbd = length(bound)
  nullcol = spzeros(nvtx, nbd)
  nullrow = spzeros(nbd, nvtx)
  Ax = a[1:nvtx, 1:nvtx]
  fx = f[1:nvtx]

  # set boundary condition
  xbd = xy[bound, 1]
  ybd = xy[bound, 2]
  if domain == 1
    bc = zero_bc(xbd, ybd)
  elseif domain == 2
    bc = obstacle_bc(xbd, ybd)
  elseif domain == 3
    bc = poiseuille_bc(xbd, ybd)
  end

  fx -= Ax[:, bound] * bc
  dA = zeros(nvtx, 1)
  dA[bound] = ones(nbd, 1)
  Ax[:, bound] = nullcol
  Ax[bound, :] = nullrow
  Ax += spdiagm(vec(dA), 0, nvtx, nvtx)
  fx[bound] = bc
  (Ax, fx)

end
