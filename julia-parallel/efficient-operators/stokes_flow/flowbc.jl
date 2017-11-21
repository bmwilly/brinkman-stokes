include("regcavity_flow.jl")
include("poiseuille_flow.jl")

###FLOWBC imposes inflow boundary conditions
# input
#   A         vector diffusion matrix
#   B         divergence matrix
#   f         velocity RHS vector
#   g         pressure RHS vector
#   xy        vertex coordinate vector
#   bound     boundary vertex vector
# output
#   Ast       vector diffusion matrix
#   Bst       divergence matrix
#   fst       velocity RHS vector
#   gst       pressure RHS vector
function flowbc(a, b, f, g, xy, bound, domain)

  nu = length(f); np = length(g)
  nvtx = Int(nu / 2); nbd = length(bound)
  nullcol = spzeros(nvtx, nbd)
  nullpcol = spzeros(np, nbd)
  Ax = a[1:nvtx, 1:nvtx]
  Ay = a[nvtx + 1:nu, nvtx + 1:nu]
  Bx = b[1:np, 1:nvtx]
  By = b[1:np, nvtx + 1:nu]
  fx = f[1:nvtx]
  fy = f[nvtx + 1:nu]
  gz = g

  # set boundary condition
  xbd = xy[bound, 1]; ybd = xy[bound, 2]
  if domain == 1
    (bcx, bcy) = regcavity_flow(xbd, ybd)
  elseif domain == 2
    (bcx, bcy) = poiseuille_flow(xbd, ybd)
  end

  # impose boundary condition
  fx -= Ax[:, bound] * bcx
  fy -= Ay[:, bound] * bcy
  gz -= Bx[:, bound] * bcx
  gz -= By[:, bound] * bcy
  dA = zeros(nvtx, 1)
  dA[bound] = ones(nbd, 1)

  Axt = Ax'
  Axt[:, bound] = nullcol
  Ax = Axt'
  Ax[:, bound] = nullcol
  Ax += spdiagm(vec(dA), 0, nvtx, nvtx)
  fx[bound] = bcx

  Ayt = Ay'
  Ayt[:, bound] = nullcol
  Ay = Ayt'
  Ay[:, bound] = nullcol
  Ay += spdiagm(vec(dA), 0, nvtx, nvtx)
  fy[bound] = bcy

  Bx[:, bound] = nullpcol
  By[:, bound] = nullpcol

  az = [Ax spzeros(nvtx, nvtx); spzeros(nvtx, nvtx) Ay]
  bz = [Bx By]
  fz = [fx; fy]

  (az, bz, fz, gz)

end
