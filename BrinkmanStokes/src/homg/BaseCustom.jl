using Base.Cartesian

@generated function indunique(A::AbstractArray{T, N}, dim::Int) where {T, N}
    return quote
        1 <= dim <= $N || return copy(A)
        hashes = zeros(Uint, size(A, dim))

        # Compute hash for each row
        k = 0
        @nloops $N i A d -> (
            if d == dim
                k = i_d
            end
        ) begin
            @inbounds hashes[k] = hash(hashes[k], hash((@nref $N A i)))
        end

        # Collect index of first row for each hash
        uniquerow = Array(Int, size(A, dim))
        firstrow = Dict{Base.Prehashed, Int}()
        for k in 1:size(A, dim)
            uniquerow[k] = get!(firstrow, Base.Prehashed(hashes[k]), k)
        end
        uniquerows = collect(values(firstrow))

        # Check for collisions
        collided = falses(size(A, dim))
        @inbounds begin
            @nloops $N i A d -> (
                if d == dim
                    k = i_d
                    j_d = uniquerow[k]
                else
                    j_d = i_d
                end
            ) begin
                if (@nref $N A j) != (@nref $N A i)
                    collided[k] = true
                end
            end
        end

        if any(collided)
            nowcollided = BitArray(size(A, dim))
            while any(collided)
                # Collect index of first row for each collided hash
                empty!(firstrow)
                for j in 1:size(A, dim)
                    collided[j] || continue
                    uniquerow[j] = get!(firstrow, Base.Prehashed(hashes[j]), j)
                end
                for v in values(firstrow)
                    push!(uniquerows, v)
                end

                # Check for collisions
                fill!(nowcollided, false)
                @nloops $N i A d -> begin
                    if d == dim
                        k = i_d
                        j_d = uniquerow[k]
                        (!collided[k] || j_d == k) && continue
                    else
                        j_d = i_d
                    end
                end begin
                    if (@nref $N A j) != (@nref $N A i)
                        nowcollided[k] = true
                    end
                end
                (collided, nowcollided) = (nowcollided, collided)
            end
        end

        sort!(uniquerows)
    end
end
