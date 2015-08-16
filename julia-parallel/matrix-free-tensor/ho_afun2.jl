push!(LOAD_PATH, "/Users/bwilliams/Documents/brinkman-stokes/julia-parallel/julia-homg")
# push!(LOAD_PATH, "/home/bmw313/Documents/brinkman-stokes/julia-parallel/julia-homg")

@everywhere using Basis
@everywhere using Mesh
@everywhere using Xform
@everywhere using Grids
@everywhere using Tensor
@everywhere using HexMeshGrids
@everywhere require("Refel.jl")

@everywhere include("hos_homg.jl")

@everywhere order = 4
@everywhere msize = 5
@everywhere dim = 2

@everywhere params = hos_homg(order, msize, dim)

@everywhere mesh = params["mesh"];
@everywhere order = params["order"];
@everywhere dof = params["dof"];
@everywhere ne = params["ne"];
@everywhere NP = params["NP"];
@everywhere bdy = params["bdy"];
@everywhere eMat = params["eMat"];

# u = SharedArray(Float64, 2dof)
# u += rand(2dof)
# w = SharedArray(Float64, length(u))
@everywhere u = rand(2dof)
@everywhere w = zeros(2dof)

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
# nchunks = length(procs())
# splits = [iround(s) for s in linspace(0, ne, nchunks + 1)]
# idx = 3
# erange = splits[idx]+1:splits[idx+1]

# @everywhere function loop_u!(mesh::Mesh.Hexmesh, order::Int, Ux::SharedArray, u::SharedArray, erange::UnitRange)
#   @show erange
#   for e in erange
#     idx = Mesh.get_node_indices(mesh, e, order)
#     Ux[:, e] = u[idx]
#   end
#   Ux
# end
#
# @everywhere function loop_u(mesh::Mesh.Hexmesh, order::Int, Ux::SharedArray, u::SharedArray, erange::UnitRange)
#   @show erange
#   for e in erange
#     idx = Mesh.get_node_indices(mesh, e, order)
#     Ux[:, e] = u[idx]
#   end
#   Ux
# end

Ux2 = zeros(NP, ne)
for e = 1:ne
  idx = Mesh.get_node_indices(mesh, e, order)
  Ux2[:, e] = u[idx]
end


@everywhere function loop_u(mesh::Mesh.Hexmesh, order::Int, dof::Int, u, Ux, Uy)
  for e = 1:size(Ux,2)
    idx = Mesh.get_node_indices(mesh, e, order)
    Ux[:, e] = u[idx]
    Uy[:, e] = u[idx+dof]
  end
  Ux,Uy
end

@everywhere loop_u_chunk(Ux, Uy) = loop_u(mesh, order, dof, u, Ux, Uy)
# @everywhere loop_u_shared!(Ux::SharedArray) = loop_u!(mesh, order, Ux, u, erange)

Ux = dzeros(NP, ne)
Uy = dzeros(NP, ne)
Ux_chunks = map(fetch, { (@spawnat p localpart(Ux)) for p = procs(Ux) })
Uy_chunks = map(fetch, { (@spawnat p localpart(Uy)) for p = procs(Uy) })
U_mats = map(loop_u_chunk, Ux_chunks, Uy_chunks)
U_matsp = pmap(loop_u_chunk, Ux_chunks, Uy_chunks)
Ux_mats = [ x[1] for x in U_mats ]
Uy_mats = [ x[2] for x in U_mats ]
Ux = reduce(hcat, Ux_mats)
Uy = reduce(hcat, Uy_mats)


# Ux,Uy = reduce(hcat, pmap(loop_u_chunk, Ux_chunks, Uy_chunks))

# p = 3
# loop_u!(mesh, order, Ux, u, erange)
# remotecall(p, loop_u, mesh, order, Ux, u, erange)
# # pmap((mesh, order, Ux, u, erange) -> loop_u!(mesh, order, Ux, u, erange), (mesh, order, Ux, u, erange))
#
# # # loop over elements
# @sync @parallel for e = 1:ne
# for e = 1:ne
#   idx = Mesh.get_node_indices(mesh, e, order)
#   Ux[:, e] = u[idx]
#   Uy[:, e] = u[idx+dof]
# end
# #
# @sync begin
#   for p in procs()
#     @async remotecall_wait(p, loop_u!, mesh, order, Ux, u, erange)
#   end
# end
#
Wx = eMat * Ux; Wy = eMat * Uy
# # @sync @parallel for e = 1:ne
# for e = 1:ne
#   idx = Mesh.get_node_indices(mesh, e, order)
#   w[idx] += Wx[:, e]
#   w[idx+dof] += Wy[:, e]
# end
#
# w[bdy] = uu[bdy]
# w[bdy+dof] = uu[bdy+dof]
