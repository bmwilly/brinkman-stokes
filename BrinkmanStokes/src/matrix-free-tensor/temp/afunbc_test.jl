include("../stokes_flow/afunbc.jl")
include("afunbc2.jl")

kparams = square_stokes()

xy = kparams["xy"]; xyp = kparams["xyp"]
nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])

u = linspace(1, nu + np, nu + np); w1 = afunbc(u, kparams)
u = linspace(1, nu + np, nu + np); w2 = afunbc3(u, kparams)
findall(abs(w1 - w2) .> 1e-5)

# tic(); for cnt = 1:100; u = rand(nu+np); w = afunbc(u, kparams); end
# etoc = toc()
#
# tic(); for cnt = 1:100; u = rand(nu+np); w = afunbc(u, kparams); end
# etoc = toc()
