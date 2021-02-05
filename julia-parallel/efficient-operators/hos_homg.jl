push!(LOAD_PATH, "$(homedir())/Documents/brinkman-stokes/julia-parallel/julia-homg")
@everywhere using ParallelSparseMatMul
@everywhere include("helpers/helper_functions.jl")
@everywhere using LinearOperators
@everywhere using ParallelSparseMatMul
@everywhere using Mesh
@everywhere using Xform
@everywhere include("../julia-homg/Refel.jl")

function hos_homg(order, msize, dim)
    nelems = [2^msize]
    m = Mesh.Hexmesh(tuple(repeat(nelems, 1, dim)...), Xform.identity)
    dof = prod([m.nelems...] * order + 1)
    K, M = Mesh.assemble_poisson(m, order)
    k1, k2 = size(K)
    A = [K spzeros(k1, k2); spzeros(k1, k2) K]
    S = share(A)
    L = operator(A)
    tic()
    for cnt = 1:100
        u = Base.shmem_rand(2dof);
        w = S * u;
    end
    etoc = toc()
end
