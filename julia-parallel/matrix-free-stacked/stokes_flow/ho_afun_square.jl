###HO_AFUN matrix-free A block operator for high order stokes
function ho_afun_square(u, params)
  # mesh = params["mesh"];
  # order = params["order"];
  dof = params["dof"];
  ne = params["ne"];
  bdy = params["bdy"];
  # refel = params["refel"];
  # mv = share(params["mv"])
  idxs = share(params["idxs"])
  eMat = share(params["eMat"])
  w = SharedArray(Float64, length(u))

  # zero dirichlet bdy conditions
  uu = copy(u)
  u[bdy] = zeros(length(bdy))
  u[bdy+dof] = zeros(length(bdy))

  # loop over elements
  # for e = 1:ne
  #   idx = Mesh.get_node_indices(mesh, e, order)
  #   pts = Mesh.element_nodes(mesh, e, refel)
  #   (detJac, Jac) = Mesh.geometric_factors(mesh, refel, pts)
  #   eMat = Mesh.element_stiffness(mesh, e, refel, detJac, Jac)
  #   w[idx] += eMat * u[idx]
  #   w[idx+dof] += eMat * u[idx+dof]
  # end

  @sync begin
    for p in procs()
      @async remotecall_wait(p, ho_afun_loop_chunk!, w, u, dof, idxs, eMat)
    end
  end

  w[bdy] = uu[bdy]
  w[bdy+dof] = uu[bdy+dof]
  vec(w)
end

@everywhere function ho_afun_loop!(w::SharedArray, u::SharedArray, dof, idxs, eMat, prange)
  pnel = length(prange)
  idxse = idxs[prange,:]
  for e = 1:pnel
    idx = vec(idxse[e,:])
    w[idx] += eMat * u[idx]
    w[idx+dof] += eMat * u[idx+dof]
  end
  w
end

@everywhere ho_afun_loop_chunk!(w, u, dof, idxs, eMat) = ho_afun_loop!(w, u, dof, idxs, eMat, myrange(idxs))

# @everywhere function ho_afun_loop!(w::SharedArray, u::SharedArray, mesh, order, dof, refel, idxs, prange)
#   pnel = length(prange)
#   # mve = mv[prange,:]
#   idxse = idxs[prange,:]
#   for e = 1:pnel
#     # idx = vec(mve[e,:]')
#     idx = vec(idxse[e,:])
#     pts = Mesh.element_nodes(mesh, e, refel)
#     detJac,Jac = Mesh.geometric_factors(mesh, refel, pts)
#     eMat = Mesh.element_stiffness(mesh, e, refel, detJac, Jac)
#     w[idx] += eMat * u[idx]
#     w[idx+dof] += eMat * u[idx+dof]
#   end
#   w
# end
#
# @everywhere ho_afun_loop_chunk!(w, u, mesh, order, dof, refel, idxs) = ho_afun_loop!(w, u, mesh, order, dof, refel, idxs, myrange(idxs))
