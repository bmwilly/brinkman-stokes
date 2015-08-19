include("hos_homg.jl")
include("mv_fun.jl")

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

@everywhere function loop_elem!(w, u, Ux, Uy, aes, mv, nvtx, erange)
  # @show erange
  for e in erange
    ind = vec(mv[e, :]')
    Ux[:, e] = u[ind]
    Uy[:, e] = u[ind+nvtx]
  end
  Wx = aes * Ux
  Wy = aes * Uy
  for e in erange
    ind = vec(mv[e, :]')
    w[ind] += Wx[:, e]
    w[ind+nvtx] += Wy[:, e]
  end
  w
end

@everywhere loop_elem_chunk!(w, u, Ux, Uy, aes, mv, nvtx) = loop_elem!(w, u, Ux, Uy, aes, mv, nvtx, myrange(mv))

msize = 9
kparams = square_stokes(msize)

xy = kparams["xy"]; xyp = kparams["xyp"]
nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
u = linspace(1, nu+np, nu+np)
xy = kparams["xy"];
xyp = kparams["xyp"];
mv = share(kparams["mv"]);
bound = kparams["bound"];
ae = kparams["ae"]

# get variables
nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
nel = length(mv[:, 1])
aes = squeeze(ae[1, :, :], 1)

# zero dirichlet boundary conditions
uu = copy(u)
u[bound] = zeros(length(bound))
u[bound+nvtx] = zeros(length(bound))

w = SharedArray(Float64, nu+np)
Ux = SharedArray(Float64, (size(mv,2), size(mv,1)))
Uy = SharedArray(Float64, (size(mv,2), size(mv,1)))

tic()
@sync begin
  for p in procs()
    @async remotecall_wait(p, loop_elem_chunk!, w, u, Ux, Uy, aes, mv, nvtx)
  end
end
etoc = toc()

# tic()


w[bound] = uu[bound]
w[bound+nvtx] = uu[bound+nvtx]
vec(w)
@show etoc
