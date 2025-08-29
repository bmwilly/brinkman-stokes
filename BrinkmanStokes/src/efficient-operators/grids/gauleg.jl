###GAULEG 1D Gauss-Legendre quadrature rule
# input
#   a,b           interval endpoints
#   n             n-point quadrature rule
# output
#   w             n-by-1 weights of quadrature rule
#   x             n-by-1 abscissae of the quadrature rule
function gauleg(a, b, n)

    # computes quadrature rule for reference interval [-1,1]
    if n == 1
        x = 0
        w = 2
    else
        c = zeros(n - 1, 1)
        for i in 1:(n - 1)
            c[i] = i / sqrt(4 * i * i - 1)
        end
        c = vec(c)
        J = diagm(c, -1) + diagm(c, 1)
        x, ev = eig(J)
        w = (2 * (ev[1, :] .* ev[1, :]))'
    end

    # if interval != [-1,1] then remap the interval and rescale the weights
    xm = (b + a) / 2 # midpoint
    xl = (b - a) / 2 # area of the requested interval

    x = xm + xl * x
    w = w * xl

    return quadrule = {"w" => w, "x" => x}

end
