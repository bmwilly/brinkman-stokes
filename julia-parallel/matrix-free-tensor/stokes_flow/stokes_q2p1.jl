# include("../diffusion/deriv.jl")
# include("../diffusion/qderiv.jl")
# include("../diffusion/lderiv.jl")
# require("diffusion/deriv.jl")
# require("diffusion/qderiv.jl")
# require("diffusion/lderiv.jl")
# include("diffusion/deriv.jl")
# include("diffusion/qderiv.jl")
# include("diffusion/lderiv.jl")
include("../diffusion/deriv.jl")
include("../diffusion/qderiv.jl")
include("../diffusion/lderiv.jl")

###STOKES_Q2P1 Q2-P1 matrix generator
#input
#    xy          Q2 nodal coordinate vector
#    xyp         Q1 nodal coordinate vector
#    mv          Q2 element mapping matrix
#output
#    ae          local Q2 diffusion derivative matrix
#    bxe         local Q2-Q1 divergence matrix
#    bye         local Q2-Q1 divergence matrix
#    ge          local Q2 vector mass matrix
#    qe          local Q1 mass matrix
#    bbxe        local Q2 x-derivative matrix
#    bbye        local Q2 y-derivative matrix
function stokes_q2p1(grid)

    xy = grid["xy"]; xyp = grid["xyp"]; mv = grid["mv"]
    nngpt = 9
    x = xy[:, 1]; y = xy[:, 2]
    xp = xyp[:, 1]; yp = xyp[:, 2]
    nvtx = length(x); nu = 2nvtx; np = 3length(xp)
    nel = length(mv[:, 1])
    mp = [[1:3:3nel]' [2:3:3nel]' [3:3:3nel]']

    println("setting up Q2-P1 matrices... ")

    ## Gauss point integration rules
    # 3x3 Gauss points
    s = zeros(nngpt, 1)
    t = zeros(nngpt, 1)
    wt = zeros(nngpt, 1)
    gpt = sqrt(0.6);
    s[1] = -gpt; t[1] = -gpt; wt[1]=25/81;
    s[2] =  gpt; t[2] = -gpt; wt[2]=25/81;
    s[3] =  gpt; t[3] =  gpt; wt[3]=25/81;
    s[4] = -gpt; t[4] =  gpt; wt[4]=25/81;
    s[5] =  0.0; t[5] = -gpt; wt[5]=40/81;
    s[6] =  gpt; t[6] =  0.0; wt[6]=40/81;
    s[7] =  0.0; t[7] =  gpt; wt[7]=40/81;
    s[8] = -gpt; t[8] =  0.0; wt[8]=40/81;
    s[9] =  0.0; t[9] =  0.0; wt[9]=64/81;

    # inner loop over elements
    xlv = zeros(nel, 4)
    ylv = zeros(nel, 4)
    for ivtx = 1:4
        xlv[:, ivtx] = x[mv[:, ivtx]]
        ylv[:, ivtx] = y[mv[:, ivtx]]
    end

    # initialize derivative matrices
    ae = zeros(nel, 9, 9)
    bxe = zeros(nel, 3, 9)
    bye = zeros(nel, 3, 9)
    ge = zeros(nel, 9, 9)
    qe = zeros(nel, 3, 3)
    bbxe = zeros(nel, 9, 9)
    bbye = zeros(nel, 9, 9)

    # loop over Gauss points
    for igpt = 1:nngpt
        sigpt = s[igpt]
        tigpt = t[igpt]
        wght = wt[igpt]

        # evaluate derivatives, etc.
        (jac, invjac, phi, dphidx, dphidy) = deriv(sigpt, tigpt, xlv, ylv)
        (psi, dpsidx, dpsidy) = qderiv(sigpt, tigpt, xlv, ylv)
        (chi, dchidx, dchidy) = lderiv(sigpt, tigpt, xlv, ylv)
        for j = 1:9
            for i = 1:9
                ae[:, i, j] += wght * dpsidx[:, i] .* dpsidx[:, j] .* invjac[:]
                ae[:, i, j] += wght * dpsidy[:, i] .* dpsidy[:, j] .* invjac[:]
                ge[:, i, j] += wght * psi[:, i] .* psi[:, j] .* jac[:]
                bbxe[:, i, j] -= wght * psi[:, i] .* dpsidx[:, j]
                bbye[:, i, j] -= wght * psi[:, i] .* dpsidy[:, j]
            end
            for i = 1:3
                bxe[:, i, j] -= wght * chi[:, i] .* dpsidx[:, j]
                bye[:, i, j] -= wght * chi[:, i] .* dpsidy[:, j]
            end
        end

        for j = 1:3
          for i = 1:3
            qe[:, i, j] += wght * chi[:, i] .* chi[:, j] .* jac[:]
          end
        end

    end # end of Gauss point loop

    println("done")

    ae = squeeze(ae[1, :, :], 1)
    bxe = squeeze(bxe[1, :, :], 1)
    bye = squeeze(bye[1, :, :], 1)
    ge = squeeze(ge[1, :, :], 1)
    qe = squeeze(qe[1, :, :], 1)
    bbxe = squeeze(bbxe[1, :, :], 1)
    bbye = squeeze(bbye[1, :, :], 1)

    elem_mats = {
      "ae" => ae,
      "bxe" => bxe,
      "bye" => bye,
      "ge" => ge,
      "qe" => qe,
      "bbxe" => bbxe,
      "bbye" => bbye
    }

end
