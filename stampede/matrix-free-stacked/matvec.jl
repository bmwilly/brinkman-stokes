using LinearOperators

reload("stokes_flow/square_stokes.jl")
# reload("stokes_flow/kfunbc.jl")
reload("stokes_flow/flop_count/kfunbc.jl")

@time kparams = square_stokes()

xy = kparams["xy"]; xyp = kparams["xyp"]
nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])

K = u -> kfunbc(u, kparams)
KK = LinearOperator(nu+np, Float64, K)

@time for cnt = 1:100; u = vec(rand(nu + np, 1)); w = KK*u; end

w,nflops = K(vec(rand(nu + np, 1)))
println("nflops: $(nflops)")
