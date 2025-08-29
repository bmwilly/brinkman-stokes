# SharedType = Union{SharedArray, SharedSparseMatrixCSC}

function input(prompt::String = "")
    print(prompt)
    return chomp(readline())
end

meshgrid(v::AbstractVector) = meshgrid(v, v)

function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where {T}
    m, n = length(vy), length(vx)
    vx = reshape(vx, 1, n)
    vy = reshape(vy, m, 1)
    return (repeat(vx, m, 1), repeat(vy, 1, n))
end

function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}, vz::AbstractVector{T}) where {T}
    m, n, o = length(vy), length(vx), length(vz)
    vx = reshape(vx, 1, n, 1)
    vy = reshape(vy, m, 1, 1)
    vz = reshape(vz, 1, 1, o)
    om = ones(Int, m)
    on = ones(Int, n)
    oo = ones(Int, o)
    return (vx[om, :, oo], vy[:, on, oo], vz[om, on, :])
end

function meshgrid(vx, vy)
    m, n = length(vy), length(vx)
    vx = reshape(vx, 1, n)
    vy = reshape(vy, m, 1)
    return (repeat(vx, m, 1), repeat(vy, 1, n))
end

function ismember(main_array, sub_array)
    out = Int8.(zeros(length(main_array)))
    match_index = findall(in(sub_array), main_array)
    out[match_index] .= 1
    return out
end

@everywhere function myrange(mv::SharedArray)
    ind = indexpids(mv)
    if ind == 0
        return 1:0
    end
    nchunks = length(procs(mv))
    splits = [iround(s) for s in linspace(0, size(mv, 1), nchunks + 1)]
    (splits[ind] + 1):splits[ind + 1]
end
