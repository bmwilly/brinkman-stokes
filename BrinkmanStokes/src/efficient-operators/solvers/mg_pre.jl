###MG_PRE presmoothing for GMG
# input
# 	A 				coefficient matrix
# 	xs0 			initial iterate
# 	f 				rhs
# 	ns 				number of smoothing steps
# 	Qs 				structure containing smoothing operators
# 	level 			grid level
# 	sweeps 			type of sweeping strategy used for Gauss-Seidel smoothing
# output
# 	xs 				result of presmoothing applied to x0
function mg_pre(A, xs0, f, ns, Qs, level, sweeps)

    xs = xs0
    if ns == 0
        return xs
    end

    L1 = Qs[level]["L1"];    U1 = Qs[level]["U1"]
    L2 = Qs[level]["L2"];    U2 = Qs[level]["U2"]
    L3 = Qs[level]["L3"];    U3 = Qs[level]["U3"]
    L4 = Qs[level]["L4"];    U4 = Qs[level]["U4"]

    r = f - A * xs

    xs += U1 \ (L1 \ r); r = f - A * xs
    if sweeps >= 2
        xs += U2 \ (L2 \ r); r = f - A * xs

        if sweeps >= 3
            xs += U3 \ (L3 \ r); r = f - A * xs

            if sweeps >= 4
                xs += U4 \ (L4 \ r)
            end
        end
    end

    if ns > 1
        k = 1
        while k < ns
            r = f - A * xs
            xs += U1 \ (L1 \ r); r = f - A * xs

            if sweeps >= 2
                xs += U2 \ (L2 \ r); r = f - A * xs

                if sweeps >= 3
                    xs += U3 \ (L3 \ r); r = f - A * xs

                    if sweeps >= 4
                        xs += U4 \ (L4 \ r)
                    end
                end
            end

            k += 1
        end
    end

    return xs

end
