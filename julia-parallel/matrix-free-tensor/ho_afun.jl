###HO_AFUN matrix-free A block operator for high order stokes
function ho_afun(u, params)

  # mesh = params["mesh"];
  # order = params["order"];
  # dof = params["dof"]; ne = params["ne"]; NP = params["NP"];
  # bdy = params["bdy"]; eMat = params["eMat"];
  # w = zeros(length(u))
  # w = SharedArray(Float64, length(u))
  # eMat = share(eMat)
  # w = share(zeros(length(u)))

  @everywhere mesh = params["mesh"];
  @everywhere order = params["order"];
  @everywhere dof = params["dof"];
  @everywhere ne = params["ne"];
  @everywhere NP = params["NP"];
  @everywhere bdy = params["bdy"];
  @everywhere eMat = params["eMat"];
  @everywhere w = zeros(2dof)

  # zero dirichlet bdy conditions
  uu = copy(u)
  u[bdy] = zeros(length(bdy))
  u[bdy+dof] = zeros(length(bdy))

  # Ux = zeros(NP, ne); Uy = zeros(NP, ne)
  # Ux = SharedArray(Float64, (NP, ne))
  # Uy = SharedArray(Float64, (NP, ne))
  # Ux = share(zeros(NP, ne)); Uy = share(zeros(NP, ne));
  Ux = dzeros(NP, ne)
  Uy = dzeros(NP, ne)

  # loop over elements
  # @sync @parallel for e = 1:ne
  # for e = 1:ne
  #   idx = Mesh.get_node_indices(mesh, e, order)
  #   Ux[:, e] = u[idx]
  #   Uy[:, e] = u[idx+dof]
  # end

  Ux = dzeros(NP, ne)
  Uy = dzeros(NP, ne)
  Ux_chunks = map(fetch, { (@spawnat p localpart(Ux)) for p = procs(Ux) })
  Uy_chunks = map(fetch, { (@spawnat p localpart(Uy)) for p = procs(Uy) })
  U_mats = map(loop_u_chunk, Ux_chunks, Uy_chunks)
  # U_matsp = pmap(loop_u_chunk, Ux_chunks, Uy_chunks)
  Ux_mats = [ x[1] for x in U_mats ]
  Uy_mats = [ x[2] for x in U_mats ]
  Ux = reduce(hcat, Ux_mats)
  Uy = reduce(hcat, Uy_mats)


  # @sync begin
  #   for p in procs()
  #     @async remotecall_wait(p, loop_u!, myrange(Ux), mesh, order, Ux, Uy, u, dof)
  #   end
  # end

  Wx = eMat * Ux; Wy = eMat * Uy
  # @sync @parallel for e = 1:ne
  for e = 1:ne
    idx = Mesh.get_node_indices(mesh, e, order)
    w[idx] += Wx[:, e]
    w[idx+dof] += Wy[:, e]
  end

  w[bdy] = uu[bdy]
  w[bdy+dof] = uu[bdy+dof]
  vec(w)
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


# @everywhere function loop_u!(erange, mesh, order, Ux, Uy, u, dof)
#   @show erange
#   for e in erange
#     idx = Mesh.get_node_indices(mesh, e, order)
#     Ux[:, e] = u[idx]
#     Uy[:, e] = u[idx+dof]
#   end
#   Ux,Uy
# end

# @everywhere function loop_shared_u!()

# @everywhere function myrange(Ux::SharedArray)
#   idx = indexpids(Ux)
#   if idx == 0
#     # this worker is not assigned a piece
#     return 1:0
#   end
#   nchunks = length(procs(Ux))
#   splits = [iround(s) for s in linspace(0, size(Ux, 2), nchunks + 1)]
#   splits[idx]+1:splits[idx+1]
# end
