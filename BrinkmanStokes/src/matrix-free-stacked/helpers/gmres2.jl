using KrylovMethods

function gmres2(A::LinearOperator, b::Vector, restrt::Int; tol::Real = 1.0e-2, maxIter::Int = 100, M::LinearOperator = opEye(length(b)), x::Vector = [], out::Int = 0)

    n = length(b)
    if isempty(x)
        x = zeros(n)
        r = M * b
    else
        r = M * (b - A * x)
    end

    bnrm2 = norm(b)
    if bnrm2 == 0.0
        bnrm2 = 1.0
    end

    err = norm(r) / bnrm2
    if err < tol
        flag = NaN; iter = NaN; resvec = NaN
        return x, flag, err, iter, resvec
    end

    return xst, flag, err, iter, resvec = gmres(
        A, b, restrt;
        tol = tol, maxIter = maxIter, M = M, x = x, out = out
    )

end
