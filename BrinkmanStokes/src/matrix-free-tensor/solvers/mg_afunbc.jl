reload("solvers/mg_afun.jl")

###MG_AFUNBC modified operator with zero boundary condition imposed
function mg_afunbc(u, mparams)

  bound = mparams["bound"]

  w = mg_afun(u, mparams)
  w[bound] = zeros(length(bound), 1)
  vec(w)

end
