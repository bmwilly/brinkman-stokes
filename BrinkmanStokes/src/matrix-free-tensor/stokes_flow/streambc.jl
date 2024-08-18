###STREAMBC imposes Dirichlet BC on the stream function
# input
#   f       RHS vector
#   xy      vertex coordinate vector
#   bound   boundary vertex vector
# output
#   f       RHS vector
function streambc(f, kparams)

  xy = kparams["xy"]; bound = kparams["bound"]

  nvtx = length(f)
  fx = f[1:nvtx]

  # set boundary condition
  xbd = xy[bound, 1]
  bcx = 0 * xbd
  fx[bound] = bcx
  fx

end
