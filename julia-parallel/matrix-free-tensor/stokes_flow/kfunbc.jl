# include("afunbc.jl")
include("ho_afun.jl")
include("bfunbc.jl")
include("btfunbc.jl")

###KFUNBC matrix-free stiffness operator
function kfunbc(u, kparams)
    w = ho_afun(u, kparams) + bfunbc(u, kparams) + btfunbc(u, kparams)
end
