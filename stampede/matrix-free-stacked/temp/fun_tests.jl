using LinearOperators
reload("stokes_flow/afun.jl")
reload("stokes_flow/afunbc.jl")
reload("stokes_flow/afunstreambc.jl")
reload("stokes_flow/agal.jl")
reload("stokes_flow/agal_diag.jl")
reload("stokes_flow/bbfunplot.jl")
reload("stokes_flow/bfun.jl")
reload("stokes_flow/bfunbc.jl")
reload("stokes_flow/btfun.jl")
reload("stokes_flow/btfunbc.jl")
reload("stokes_flow/gfun.jl")
reload("stokes_flow/ho_afun.jl")
reload("stokes_flow/kfunbc.jl")
reload("stokes_flow/mfun.jl")
reload("stokes_flow/qfun.jl")
reload("stokes_flow/qfun_diag.jl")

reload("stokes_flow/square_stokes.jl")

kparams = square_stokes()

xy = kparams["xy"]; xyp = kparams["xyp"]
nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])

u = linspace(1, nu+np, nu+np)

# w = kfunbc(u, kparams)
w = mfun(u, kparams)
# w = qfun(u, kparams)
# w = gfun(u, kparams)
