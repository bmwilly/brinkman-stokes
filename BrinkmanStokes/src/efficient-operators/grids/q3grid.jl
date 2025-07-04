###Q2GRID produces grid of 16 node quadrilaterals
# input
#   nelemx          number of elements along x direction
#   nelemy          number of elements along y direction
# output
#   element_node    nodes that form each element
function q3grid(nelemx, nelemy)

    element = 0
    element_node = zeros(Int, 16, nelemx * nelemy)

    for j = 1:nelemy
        for i = 1:nelemx
            base = (j - 1) * 3 * (3 * nelemx + 1) + 3i - 2
            element += 1

            element_node[1, element] = base
            element_node[2, element] = base + 1
            element_node[3, element] = base + 2
            element_node[4, element] = base + 3
            element_node[5, element] = base + (3nelemx + 1)
            element_node[6, element] = base + (3nelemx + 1) + 1
            element_node[7, element] = base + (3nelemx + 1) + 2
            element_node[8, element] = base + (3nelemx + 1) + 3
            element_node[9, element] = base + 2 * (3nelemx + 1)
            element_node[10, element] = base + 2 * (3nelemx + 1) + 1
            element_node[11, element] = base + 2 * (3nelemx + 1) + 2
            element_node[12, element] = base + 2 * (3nelemx + 1) + 3
            element_node[13, element] = base + 3 * (3nelemx + 1)
            element_node[14, element] = base + 3 * (3nelemx + 1) + 1
            element_node[15, element] = base + 3 * (3nelemx + 1) + 2
            element_node[16, element] = base + 3 * (3nelemx + 1) + 3

        end
    end

    element_node'
end
