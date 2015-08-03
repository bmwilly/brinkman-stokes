using LinearOperators
reload("stokes_flow/square_stokes.jl")
reload("stokes_flow/kfunbc.jl")
reload("stokes_flow/afunbc.jl")

function mv_fun(msize)
  kparams = square_stokes()

  xy = kparams["xy"]; xyp = kparams["xyp"]
  nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])

  tic(); for cnt = 1:100; u = rand(nu+np); w = afunbc(u, kparams); end
  etoc = toc()

  # K = u -> kfunbc(u, kparams)
  # KK = LinearOperator(nu+np, Float64, K)
  #
  # tic()
  # for cnt = 1:100; u = vec(rand(nu + np, 1)); w = KK*u; end
  # etoc = toc()
end
