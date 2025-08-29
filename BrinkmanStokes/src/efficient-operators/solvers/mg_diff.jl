include("mg_prolong.jl")
include("mg_diff_setup.jl")
include("mg_smooth.jl")

###MG_DIFF GMG preconditioner for diffusion problem
function mg_diff(x, y, Agal)

    nc = Int(log2(length(y) - 1))

    # compute new MG data
    h = 2^(1 - Float64(nc))

    println("Setting up MG data...")

    # top level
    mgdata = Array{Dict}(undef, Int(nc))
    mgdata[nc] = Dict(
        "matrix" => Agal,
        "prolong" => mg_prolong(2^nc, 2^nc, x, y)
    )

    xn = x; yn = y
    # loop over remaining levels
    for level in (nc - 1):-1:2
        xn = xn[1:2:end]; yn = yn[1:2:end]
        mgdata[level] = Dict(
            "matrix" => mg_diff_setup(xn, yn),
            "prolong" => mg_prolong(2^level, 2^level, xn, yn)
        )
    end

    println("done")

    # MG parameters
    # smooth = user_input("Smoother? (1:Jacobi, 2:Gauss-Seidel, 3:ILU)")
    smooth = 1 # Jacobi
    sweeps = 1; stype = 1

    # pre- and post-smoothing steps
    # npre = 1; npost = 1
    npre = user_input("Number of pre-smoothing steps: ")
    npost = user_input("Number of post-smoothing steps: ")

    # construct smoother
    smooth_data = mg_smooth(mgdata, nc, sweeps, smooth, stype)

    return (mgdata, smooth_data, sweeps, stype, npre, npost, nc)

end
