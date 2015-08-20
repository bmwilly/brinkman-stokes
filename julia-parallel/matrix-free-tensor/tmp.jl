using ParallelSparseMatMul
include("stokes_flow/square_stokes.jl")

@everywhere function myrange(mv::SharedArray)
  ind = indexpids(mv)
  if ind == 0
    # this worker is not assigned a piece
    return 1:0
  end
  nchunks = length(procs(mv))
  splits = [iround(s) for s in linspace(0, size(mv, 1), nchunks + 1)]
  splits[ind]+1:splits[ind+1]
end

@everywhere function loop_elem!(w::SharedArray, u::SharedArray, aes::SharedArray, mv::SharedArray, nvtx, prange::UnitRange)
  @show prange
  pnel = length(prange)
  mve = mv[prange,:]
  Ux = zeros(size(mve,2), pnel)
  Uy = zeros(size(mve,2), pnel)
  for e in 1:pnel
    ind = vec(mve[e, :]')
    Ux[:, e] = u[ind]
    Uy[:, e] = u[ind+nvtx]
  end
  Wx = aes * Ux
  Wy = aes * Uy
  for e in 1:pnel
    ind = vec(mve[e, :]')
    w[ind] += Wx[:, e]
    w[ind+nvtx] += Wy[:, e]
  end

  # for e in erange
  #   ind = vec(mv[e, :]')
  #   Ux[:, e] = u[ind]
  #   Uy[:, e] = u[ind+nvtx]
  # end
  # Wx = aes * Ux
  # Wy = aes * Uy
  # for e in erange
  #   ind = vec(mv[e, :]')
  #   w[ind] += Wx[:, e]
  #   w[ind+nvtx] += Wy[:, e]
  # end

  w
end

@everywhere loop_elem_chunk!(w, u, aes, mv, nvtx) = loop_elem!(w, u, aes, mv, nvtx, myrange(mv))

msize = 9
kparams = square_stokes(msize)
#
xy = kparams["xy"]; xyp = kparams["xyp"]
nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
xy = kparams["xy"];
xyp = kparams["xyp"];
bound = kparams["bound"];
ae = kparams["ae"]

# get variables
mv = share(kparams["mv"]);
nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
nel = length(mv[:, 1])
# aes = ae
aes = share(ae)

# zero dirichlet boundary conditions
u = share(linspace(1, nu+np, nu+np))
uu = copy(u)
u[bound] = zeros(length(bound))
u[bound+nvtx] = zeros(length(bound))

w = SharedArray(Float64, nu+np)
# Ux = SharedArray(Float64, (size(mv,2), size(mv,1)))
# Uy = SharedArray(Float64, (size(mv,2), size(mv,1)))

tic()
@sync begin
  for p in procs()
    @async remotecall_wait(p, loop_elem_chunk!, w, u, aes, mv, nvtx)
  end
end
t2 = toc()

# tic()
w[bound] = uu[bound]
w[bound+nvtx] = uu[bound+nvtx]
vec(w)

@show t2
