module TestGrid

include("Xform.jl")
include("HexMesh.jl")
include("Grids.jl")
include("HexMeshGrids.jl")

# mu = (x,y)->(1 + 1e6*( cos(2*pi*x)^2 + cos(2*pi*y)^2 ))
#
# g = HexMeshGrids.create_hexmesh_grids(2, mu, Xform.shell, [1 2 4], [8 16 32]);
#
# # (u, rr, iter) = Grids.solve(g, 150, "jacobi", 3,3, g.L, Grids.get_u0(g));
# (u, rr, iter) = Grids.solve_pcg(g, 150, "jacobi", 3,3, g.L, Grids.get_u0(g));

mu = (x,y,z) -> (1 + 1e6*( cos(2pi*x)^2 + cos(2pi*y)^2 + cos(2pi*z)^2 ))
g = HexMeshGrids.create_hexmesh_grids(3, mu, Xform.shell, [1 2 4], [2 4 8]);

# (u, rr, iter) = Grids.solve(g, 150, "jacobi", 3,3, g.L, Grids.get_u0(g));
(u, rr, iter) = Grids.solve_pcg(g, 150, "jacobi", 3,3, g.L, Grids.get_u0(g));

end
