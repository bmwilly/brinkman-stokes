reload("stokes_flow/square_stokes.jl")
reload("solvers/mg_prolong.jl")
reload("solvers/mg_diff_setup.jl")

## solve_stokes_mg
(A, B, G, Q, Bx, By, f, g, x, y, xy, xyp, bound) = square_stokes()

# boundary conditions
println("imposing (enclosed flow) boundary conditions ...")
(Ast, Bst, fst, gst) = flowbc(A, B, f, g, xy, bound)
np = length(gst)

## STOKES problem

# set parameters
tol = 1e-6; 
maxit = 100;

# specify coefficient matrix
K = [Ast Bst'; Bst spzeros(np, np)]

# construct preconditioner
println("multigrid preconditioning...")

nv = size(Ast, 1); np = size(Q, 1); nu = nv/2
Agal = Ast[1:nu, 1:nu]




## mg_diff
nc = log2(length(y) - 1)

# compute new MG data
h = 2^(1 - nc)

println("Setting up MG data...")

# top level
mgdata = Array(Dict{}, int(nc))
mgdata[nc] = {
  "matrix" => Agal, 
  "prolong" => mg_prolong(2^nc, 2^nc, x, y)
}

xn = x; yn = y

level = (nc-1)
xn = xn[1:2:end]; yn = yn[1:2:end]
# matrix = mg_diff_setup(xn, yn)




## mg_diff_setup
x = xn; y = yn;
n = length(x) - 1; np = n/2; nq = n/4
  println("np $(np)")
  nel = int(np^2)
nvtx = (n + 1) * (n + 1)
(X, Y) = meshgrid(x, y)
xx = reshape(X', nvtx, 1)
yy = reshape(Y', nvtx, 1)
xy = [xx[:] yy[:]]

# assemply process
kx = 1; ky = 1; mel = 0
mv = zeros(Int64, nel, 9)
  mp = zeros(Int64, nel, 4)
  for j = 1:np
      for i = 1:np
          mref = (n + 1) * (ky - 1) + kx
          pref = (np + 1) * (j - 1) + i
          mel += 1
          nvv = zeros(9)
          nvv[1] = mref
          nvv[2] = mref + 2
          nvv[3] = mref + 2n + 4
          nvv[4] = mref + 2n + 2
          nvv[5] = mref + 1
          nvv[6] = mref + n + 3
          nvv[7] = mref + 2n + 3
          nvv[8] = mref + n + 1
          nvv[9] = mref + n + 2
          npp = zeros(4)
          npp[1] = pref
          npp[2] = pref + 1
          npp[3] = pref + np + 2
          npp[4] = pref + np + 1
          mv[mel, 1:9] = nvv[1:9]
          mp[mel, 1:4] = npp[1:4]
          kx += 2
      end
      ky += 2
      kx = 1
  end

  # compute boundary vertices and edges
  # four boundary edges
  k1 = find(xy[:,2] .== -1)
  e1 = []
  for k = 1:mel
      if any(mv[k,5] .== k1)
          e1 = [e1, k]
      end
  end
  ef1 = ones(size(e1))

  k2 = find((xy[:,1] .== 1) & (xy[:,2] .< 1) & (xy[:,2] .> -1))
  e2 = []
  for k = 1:mel
      if any(mv[k,6] .== k2)
          e2 = [e2, k]
      end
  end
  ef2 = 2*ones(size(e2))

  k3 = find(xy[:,2] .== 1)
  e3 = []
  for k = 1:mel
      if any(mv[k,7] .== k3)
          e3 = [e3, k]
      end
  end
  ef3 = 3*ones(size(e3))

  k4 = find((xy[:,1] .== -1) & (xy[:,2] .< 1) & (xy[:,2] .> -1))
  e4 = []
  for k = 1:mel
      if any(mv[k,8] .== k4)
          e4 = [e4, k]
      end
  end
  ef4 = 4*ones(size(e4))

  bound = sort([k1; k2; k3; k4])
  mbound = [e1' ef1'; e2' ef2'; e3' ef3'; e4' ef4']

  # set up matrices for Q1 approximation
  println("mv $(size(mv))")
  # (ev, ebound) = mg_q1grid(x, y, xy, mv, bound, mbound)





## mg_q1grid
xx = xy[:, 1]; yy = xy[:, 2]; nvtx = length(xx)
adj = spzeros(nvtx, nvtx)
mel = length(mv[:, 1]); nel = 4mel
ev = zeros(nel, 4)
ebound = zeros(16, 2)

# loop over macroelements
for k = 1:mel

  # first element
  ke = 4k - 3
  ev[ke, 1] = mv[k, 1]
  ev[ke, 2] = mv[k, 5]
  ev[ke, 3] = mv[k, 9]
  ev[ke, 4] = mv[k, 8]

  # second element
  ke = 4k - 2
  ev[ke, 1] = mv[k, 5]
  ev[ke, 2] = mv[k, 2]
  ev[ke, 3] = mv[k, 6]
  ev[ke, 4] = mv[k, 9]

  # third element
  ke = 4k - 1
  ev[ke, 1] = mv[k, 9]
  ev[ke, 2] = mv[k, 6]
  ev[ke, 3] = mv[k, 3]
  ev[ke, 4] = mv[k, 7]

  # fourth element
  ke = 4k
  ev[ke, 1] = mv[k, 8]
  ev[ke, 2] = mv[k, 9]
  ev[ke, 3] = mv[k, 7]
  ev[ke, 4] = mv[k, 4]

end

# define element edges
ect = 1
# bottom boundary edges
k1 = find(mbound[:, 2] .== 1)'
for k = mbound[k1]
  ebound[ect, 1] = 4k - 3
  ebound[ect + 1, 1] = 4k - 2
  ebound[ect, 2] = 1
  ebound[ect + 1, 2] = 1
  ect += 2
end

# right boundary edges
k2 = find(mbound[:, 2] .== 2)'
for k = mbound[k2]
  ebound[ect, 1] = 4k - 2
  ebound[ect + 1, 1] = 4k - 1
  ebound[ect, 2] = 2
  ebound[ect + 1, 2] = 2
  ect += 2
end

# top boundary edges
k3 = find(mbound[:, 2] .== 3)'
for k = mbound[k3]
  ebound[ect, 1] = 4k - 1
  ebound[ect + 1, 1] = 4k
  ebound[ect, 2] = 3
  ebound[ect + 1, 2] = 3
  ect += 2
end

# left boundary edges
k4 = find(mbound[:, 2] .== 4)'
for k = mbound[k4]
  ebound[ect, 1] = 4k
  ebound[ect + 1, 1] = 4k - 3
  ebound[ect, 2] = 4
  ebound[ect + 1, 2] = 4
  ect += 2
end

(ev, ebound)





  # (A, M, fdummy) = mg_q1diff(xy, ev)





## mg_q1diff
x = xy[:, 1]; y = xy[:, 2]
nvtx = length(x)
nel = length(ev[:, 1])
lx = maximum(x) - minimum(x); ly = maximum(y) - minimum(y)
hx = maximum(diff(x)); hy = maximum(diff(y))

println(size(x))
println(size(ev))
# println(x[ev[:,1]])
# println(ev)

# initialize global matrices
a = spzeros(nvtx, nvtx)
r = spzeros(nvtx, nvtx)
f = zeros(nvtx, 1)
ae = zeros(nel, 4, 4)

# set up 2x2 Gauss points
s = zeros(4, 1); t = zeros(4, 1)
gpt = 1.0/sqrt(3.0)
s[1] = -gpt;  t[1] = -gpt;
s[2] =  gpt;  t[2] = -gpt;
s[3] =  gpt;  t[3] =  gpt;
s[4] = -gpt;  t[4] =  gpt;  

    # inner loop over elements
    xlv = zeros(nel, 4); ylv = zeros(nel, 4)
    for ivtx = 1:4
        xlv[:, ivtx] = x[ev[:, ivtx]]
        ylv[:, ivtx] = y[ev[:, ivtx]]
    end

    # loop over 2x2 Gauss points
    for igpt = 1:4
      sigpt = s[igpt]
      tigpt = t[igpt]

      # evaluate derivatives, etc
      (jac, invjac, phi, dphidx, dphidy) = deriv(sigpt, tigpt, xlv, ylv)
      for j = 1:4
        for i = 1:4
          ae[:, i, j] += dphidx[:, i] .* dphidx[:, j] .* invjac[:]
          ae[:, i, j] += dphidy[:, i] .* dphidy[:, j] .* invjac[:]
        end
      end
    end # Gauss point loop 

    # assemble global matrix and source vector
    for krow = 1:4
      nrow = ev[:, krow]
      for kcol = 1:4
        ncol = ev[:, kcol]
        a += sparse(nrow, ncol, ae[:, krow, kcol], nvtx, nvtx)
      end
    end



















