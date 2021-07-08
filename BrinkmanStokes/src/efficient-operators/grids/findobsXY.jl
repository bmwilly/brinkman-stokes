###FINDOBSXY finds coords of grid points inside obstacle
# input
#   obs         indices in bndxy of nodes determining object boundary
#   X           array of horizontal grid coordinates
#   Y           array of vertical grid coordinates
#   bndxy       coordinates of nodes defining domain boundary
# output
#   KK          indices of mesh points interior to object
function findobsXY(obs, X, Y, bndxy)

    KK = Int[]
    if size(obs, 1) > 0
        for i = 1:size(obs, 1)
            xl = minimum(bndxy[:, 1][obs[i, :]])
            xr = maximum(bndxy[:, 1][obs[i, :]])
            yb = minimum(bndxy[:, 2][obs[i, :]])
            yt = maximum(bndxy[:, 2][obs[i, :]])
            kk = findall((X .> xl) .& (X .< xr) .& (Y .< yt) .& (Y .> yb))
            KK = vcat(KK, kk)
        end
    end

    KK
end
