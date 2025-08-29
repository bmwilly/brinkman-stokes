# require("stokes_flow/afunbc.jl")
# require("stokes_flow/bfunbc.jl")
# require("stokes_flow/btfunbc.jl")
# reload("stokes_flow/afunbc.jl")
# reload("stokes_flow/bfunbc.jl")
# reload("stokes_flow/btfunbc.jl")
reload("stokes_flow/flop_count/afunbc.jl")
reload("stokes_flow/flop_count/bfunbc.jl")
reload("stokes_flow/flop_count/btfunbc.jl")

###KFUNBC matrix-free stiffness operator
function kfunbc(u, kparams)

    xy = kparams["xy"]; xyp = kparams["xyp"]
    nvtx = length(xy[:, 1]); nu = 2nvtx; np = 3length(xyp[:, 1])

    A = x -> afunbc(x, kparams)
    B = x -> bfunbc(x, kparams)
    Bt = x -> btfunbc(x, kparams)

    # A = LinearOperator(nu+np, Float64, A)
    # B = LinearOperator(nu+np, Float64, B)
    # Bt = LinearOperator(nu+np, Float64, Bt)

    wa, na = A(u); wb, nb = B(u); wbt, nbt = Bt(u)
    w = wa + wb + wbt
    nflops = na + nb + nbt
    return w, nflops

    # wa = afunbc(u, kparams);
    # wb = bfunbc(u, kparams);
    # wbt = btfunbc(u, kparams)
    # w = wa + wb + wbt

    # wa = @spawn afunbc(u, xy, xyp, mv, bound, ae)
    # wb = @spawn bfunbc(u, xy, xyp, mv, bound, bxe, bye)
    # wbt = @spawn btfunbc(u, xy, xyp, mv, bound, bxe, bye)
    # w = fetch(wa) + fetch(wb) + fetch(wbt)

    # tic()
    # wa = afunbc(u, xy, xyp, mv, bound, ae)
    # wa = @spawn afunbc(u, xy, xyp, mv, bound, ae)
    # println("afunbc: $(toc())")

    # tic()
    # wb = bfunbc(u, xy, xyp, mv, bound, bxe, bye)
    # wb = @spawn bfunbc(u, xy, xyp, mv, bound, bxe, bye)
    # println("bfunbc: $(toc())")

    # tic()
    # wbt = btfunbc(u, xy, xyp, mv, bound, bxe, bye)
    # wbt = @spawn btfunbc(u, xy, xyp, mv, bound, bxe, bye)
    # println("btfunbc: $(toc())")

    # tic()
    # wa = convert(Array, wa)
    # wb = convert(Array, wb)
    # wbt = convert(Array, wbt)
    # w = wa + wb + wbt
    # w = fetch(wa) + fetch(wb) + fetch(wbt)
    # w = vec(w)
    # w = vec(convert(Array, w))
    # println("convert: $(toc())")

    # wa = @sync afunbc(u, xy, xyp, mv, bound, ae)
    # wb = @sync bfunbc(u, xy, xyp, mv, bound, bxe, bye)
    # wbt = @sync btfunbc(u, xy, xyp, mv, bound, bxe, bye)
    # w = wa + wb + wbt

    # wa = afunbc(u, xy, xyp, mv, bound, ae)
    # wb = bfunbc(u, xy, xyp, mv, bound, bxe, bye)
    # wbt = btfunbc(u, xy, xyp, mv, bound, bxe, bye)
    # # w = wa + wb + wbt
    # w = fetch(wa) + fetch(wb) + fetch(wbt)

end
