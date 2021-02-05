###FINDBOUNDARY generates vector bound and matrix bound
# input
#   bndxy       coordinates of nodes defining domain boundary
#   bnde        indices in bndxy of boundary edges
#   xy          vertex coordinates
#   mv          q2 element mapping matrix
# output
#   bound       dirichlet boundary vertices
#   mbound      dirichlet boundary elements
function findboundary(bndxy, bnde, xy, mv)

    bound = Int[]; mbound = [];
    for i = 1:size(bnde, 1)
        if bnde[i, 3] == 1
            if bndxy[bnde[i, 1], 1] .== bndxy[bnde[i, 2], 1]
                yl = min(bndxy[bnde[i, 1], 2], bndxy[bnde[i, 2], 2])
                yu = max(bndxy[bnde[i, 1], 2], bndxy[bnde[i, 2], 2])
                k = findall((xy[:, 1] .== bndxy[bnde[i, 1], 1]) .& (xy[:, 2] .< yu) .& (xy[:, 2] .> yl))
                bound = vcat(bound, k)
                for j = 1:size(mv, 1)
                    if any(mv[j, 6] .== k)
                        if mbound == []
                            mbound = [j 2]
                        else
                            mbound = [mbound; j 2]
                        end
                    end
                    if any(mv[j, 8] .== k)
                        if mbound == []
                            mbound = [j 4]
                        else
                            mbound = [mbound; j 4]
                        end
                    end
                end
            else
                xl = min(bndxy[bnde[i, 1], 1], bndxy[bnde[i, 2], 1])
                xu = max(bndxy[bnde[i, 1], 1], bndxy[bnde[i, 2], 1])
                k = findall((xy[:, 2] .== bndxy[bnde[i, 1], 2]) .& (xy[:, 1] .<= xu) .& (xy[:, 1] .>= xl))
                bound = vcat(bound, k)
                for j = 1:size(mv, 1)
                    if any(mv[j, 5] .== k)
                        if mbound == []
                            mbound = [j 1]
                        else
                            mbound = [mbound; j 1]
                        end
                    end
                    if any(mv[j, 7] .== k)
                        if mbound == []
                            mbound = [j 3]
                        else
                            mbound = [mbound; j 3]
                        end
                    end
                end
            end
        end
    end

    bound = sort(bound)
    (bound, mbound)

end
