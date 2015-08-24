push!(LOAD_PATH, "/Users/bwilliams/Documents/brinkman-stokes/julia-parallel/julia-homg")

@everywhere using Basis
@everywhere using Mesh
@everywhere using Xform
@everywhere using Grids
@everywhere using Tensor
@everywhere using HexMeshGrids
@everywhere require("Refel.jl")

@everywhere include("hos_homg.jl")

params = hos_homg(2, 4, 2)

mesh = params["mesh"];
order = params["order"];
dof = params["dof"];
ne = params["ne"];
NP = params["NP"];
bdy = params["bdy"];
eMat = params["eMat"];

# u = SharedArray(Float64, 2dof)
# u += rand(2dof)
# u = drand(2dof)
u = rand(2dof)
@everywhere u = rand(2dof)
# w = SharedArray(Float64, length(u))
# w = dzeros(2dof)
w = zeros(2dof)

# zero dirichlet bdy conditions
uu = copy(u)
u[bdy] = zeros(length(bdy))
u[bdy+dof] = zeros(length(bdy))

# Ux = SharedArray(Float64, (NP, ne))
# Uy = SharedArray(Float64, (NP, ne))
Ux = dzeros(NP, ne)
Uy = dzeros(NP, ne)

# @everywhere nchunks = length(procs())
# @everywhere splits = [iround(s) for s in linspace(0, ne, nchunks + 1)]
# @everywhere i = 3
# @everywhere erange = splits[i]+1:splits[i+1]
nchunks = length(procs())
splits = [iround(s) for s in linspace(0, ne, nchunks + 1)]
idx = 3
erange = splits[idx]+1:splits[idx+1]

# # loop over elements
@sync @parallel for e = 1:ne
# for e = 1:ne
  idx = Mesh.get_node_indices(mesh, e, order)
  Ux[:, e] = u[idx]
  Uy[:, e] = u[idx+dof]
end

# u_chunks = map(fetch, { (@spawnat p localpart(u)) for p = procs(u) })
Ux_chunks = map(fetch, { (@spawnat p localpart(Ux)) for p = procs(Ux) })
U = pmap(loop_u_chunk, Ux_chunks)

@sync begin
  for (i,w) in enumerate(workers())
    irange,jrange = Ux.indexes[i]
    Uxw = fetch(Ux.chunks[i])
    for e in jrange
      idx = Mesh.get_node_indices(mesh, e, order)
      Uxw[:,e] = u[idx]
    end
  end
end

#
# # @sync begin
# #   for p in procs()
# #     @async remotecall_wait(p, loop_u!, myrange(Ux), mesh, order, Ux, Uy, u, dof)
# #   end
# # end
#
# Wx = eMat * Ux; Wy = eMat * Uy
# # @sync @parallel for e = 1:ne
# for e = 1:ne
#   idx = Mesh.get_node_indices(mesh, e, order)
#   w[idx] += Wx[:, e]
#   w[idx+dof] += Wy[:, e]
# end
#
# w[bdy] = uu[bdy]
# w[bdy+dof] = uu[bdy+dof]
