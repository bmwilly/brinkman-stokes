using PyPlot
include("../hos_homg.jl")

order = int(input("Order: "))
msize = int(input("Mesh size: "))
dim = int(input("Dimension: "))

nelems = [2^msize]
m = Mesh.Hexmesh(tuple(repeat(nelems, 1, dim)...), Xform.identity)
dof = prod([m.nelems...] * order .+ 1)
K, M = Mesh.assemble_poisson(m, order)
k1, k2 = size(K)
A = [K spzeros(k1, k2); spzeros(k1, k2) K]

figure(); spy(A)

nnzs = [nnz(A[i,:]) for i = 1:size(A, 1)]
@show mean(nnzs)
@show mean(nnzs) / size(A, 1)
@show dof * 2
