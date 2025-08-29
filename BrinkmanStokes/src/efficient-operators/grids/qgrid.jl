###QGRID produces a grid of (p+1)*(p+1) node quadrilaterals
# input
#   nelemx          number of elements along x direction
#   nelemy          number of elements along y direction
#   p               order
# output
#   element_node    nodes that form each element
function qgrid(nelemx, nelemy, p)

    element = 0
    nnodes = (p + 1) * (p + 1)
    element_node = zeros(Int, nnodes, nelemx * nelemy)

    for j in 1:nelemy
        for i in 1:nelemx
            base = p * (j - 1) * (p * nelemx + 1) + p * (i - 1) + 1
            element += 1

            n = 1
            for n1 in 0:p
                for n2 in 0:p
                    element_node[n, element] = base + n1 * (p * nelemx + 1) + n2
                    n += 1
                end
            end

        end
    end

    return element_node'
end
