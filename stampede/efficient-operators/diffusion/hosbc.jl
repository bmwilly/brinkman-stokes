reload("stokes_flow/regcavity_flow.jl")
###HOSBC imposes inflow boundary condition
function hosbc(a, f, xy, bound)
  nu = length(f);
  nvtx = int(nu / 2); nbd = length(bound)
  nullcol = spzeros(nvtx, nbd)
  Ax = a[1:nvtx, 1:nvtx]
  Ay = a[nvtx + 1:nu, nvtx + 1:nu]
  fx = f[1:nvtx]
  fy = f[nvtx + 1:nu]

  # set boundary condition
  xbd = xy[bound, 1]; ybd = xy[bound, 2]
  (bcx, bcy) = regcavity_flow(xbd, ybd)

  # impose boundary condition
  fx -= Ax[:, bound] * bcx
  fy -= Ay[:, bound] * bcy
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

  az = [Ax spzeros(nvtx, nvtx); spzeros(nvtx, nvtx) Ay]
  fz = [fx; fy]

  (az, fz)
end
