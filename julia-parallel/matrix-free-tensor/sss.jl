using ParallelSparseMatMul
include("stokes_flow/solve_stokes.jl")
include("graphs/flowplot.jl")
include("helpers/helper_functions.jl")

domain = int(input("Choose problem (1/lid-driven cavity, 2/brinkman): "))
msize = int(input("Mesh size: "))
# msize = 7
# domain = 2
(sol, kparams) = solve_stokes(domain, msize)
println("done")
doplot = input("Create plot? (y/n): ")

if doplot == "y"
    flowplot(sol, kparams)
end

# x = kparams["x"]; y = kparams["y"];
# xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"]
# nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
#
# u = sol[1:nu]
# p = sol[nu+1:end]
# ux = reshape(u[1:nvtx], length(x), length(y))'
# uy = reshape(u[nvtx+1:end], length(x), length(y))'

# DelimitedFiles.writedlm("$(homedir())/Documents/brinkman-stokes/julia-parallel/matrix-free-tensor/temp/sol/x.csv", x)
# DelimitedFiles.writedlm("$(homedir())/Documents/brinkman-stokes/julia-parallel/matrix-free-tensor/temp/sol/y.csv", y)
# DelimitedFiles.writedlm("$(homedir())/Documents/brinkman-stokes/julia-parallel/matrix-free-tensor/temp/sol/ux.csv", ux)
# DelimitedFiles.writedlm("$(homedir())/Documents/brinkman-stokes/julia-parallel/matrix-free-tensor/temp/sol/uy.csv", uy)
# DelimitedFiles.writedlm("$(homedir())/Documents/brinkman-stokes/julia-parallel/matrix-free-tensor/temp/sol/kp.csv", kp)
