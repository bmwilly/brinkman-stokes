using LinearOperators
include("../../julia-homg/Basis.jl")
include("../../julia-homg/Hexmesh.jl")
include("../../julia-homg/Xform.jl")
include("../../julia-homg/Grids.jl")
include("../../julia-homg/Tensor.jl")
include("../../julia-homg/Refel.jl")
# reload("stokes_flow/ho_afun.jl")
reload("helpers/helper_functions.jl")

dim = 2
# order = int(input("Polynomial order: "))
msize = int(input("Mesh size: "))
# msize = 2
nelems = [2^msize]

m = Mesh.Hexmesh(tuple(repeat(nelems, 1, dim)...), Xform.identity)
Mesh.plot(m)
