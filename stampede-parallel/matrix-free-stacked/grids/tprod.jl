###TPROD 2D tensor-product Gauss quadrature
# input
#   quadrule          1D Gauss quadrature
# output
#   quadrule          2D Gauss quadrature
function tprod(quadrule)

  nGauss_1D = size(quadrule["w"], 1)
  nGauss_2D = nGauss_1D^2

  w = zeros(nGauss_2D, 1)
  x = zeros(nGauss_2D, 2)

  # build 2D quadrature on the reference square
  k = 1
  for i = 1:nGauss_1D
    for j = 1:nGauss_1D
      w[k] = quadrule["w"][i] * quadrule["w"][j]
      x[k, 1] = quadrule["x"][i]
      x[k, 2] = quadrule["x"][j]
      k += 1
    end
  end

  quadrule2D = {"w" => w, "x" => x}
end
