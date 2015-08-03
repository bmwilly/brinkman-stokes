reload("stokes_flow/qfun.jl")
reload("stokes_flow/gfun.jl")

function mfun(u, kparams)

  xy = kparams["xy"]; xyp = kparams["xyp"]; mv = kparams["mv"]; bound = kparams["bound"]
  nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])
  nel = length(mv[:, 1])

  wg = gfun(u, kparams); wq = qfun(u, kparams)
  w = wg + [zeros(nu); wq]
  vec(w)
end
