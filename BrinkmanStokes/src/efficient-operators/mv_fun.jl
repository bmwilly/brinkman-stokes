using LinearOperators
reload("stokes_flow/square_stokes.jl")
reload("stokes_flow/flowbc.jl")

function mv_fun(msize)

    tic()
    mats = square_stokes()
    A = mats["A"]
    B = mats["B"]
    Bx = mats["Bx"]
    By = mats["By"]
    f = mats["f"]
    g = mats["g"]
    xy = mats["xy"]
    xyp = mats["xyp"]
    bound = mats["bound"]
    x = mats["x"]
    y = mats["y"]
    Q = mats["Q"]

    # boundary conditions
    (Ast, Bst, fst, gst) = flowbc(A, B, f, g, xy, bound, 1)
    np = length(gst)
    rhs = vec([fst; gst])

    ## compute solution
    K = [Ast Bst'; Bst spzeros(np, np)]
    n, m = size(K)
    # K = LinearOperator(K)
    toc()

    tic()
    for cnt = 1:100
        u = rand(size(K, 1))
        w = K * u
    end
    etoc = toc()

    # nza = float(bool(K))
    # nzb = float(bool(vec(rand(size(K, 1), 1))))
    # f = nza * nzb
    # f = 2f - float(bool(f))
    # f = sum(sum(f))

    return (etoc)
end
