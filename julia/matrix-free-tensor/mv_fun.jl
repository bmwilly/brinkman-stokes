reload("stokes_flow/square_stokes.jl")
reload("stokes_flow/kfunbc.jl")
reload("stokes_flow/afunbc.jl")

function mv_fun(msize)
  kparams = square_stokes(msize)

  xy = kparams["xy"]; xyp = kparams["xyp"]
  nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])

  tic(); for cnt = 1:100; u = rand(nu+np); w = afunbc(u, kparams); end
  etoc = toc()
end
