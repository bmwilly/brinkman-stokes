###Q2GRID produces grid of 9 node quadrilaterals
# input
#   nelemx          number of elements along x direction
#   nelemy          number of elements along y direction
# output
#   element_node    nodes that form each element
#
#  Node labeling:
#
#    NW----N----NE
#     |          |
#     W    C     E
#     |          |
#    SW----S----SE
#
function q2grid(nelemx, nelemy)

    element = 0
    element_node = zeros(Int, 9, nelemx * nelemy)

    for j in 1:nelemy
        for i in 1:nelemx
            sw = 2 * (j - 1) * (2nelemx + 1) + 2 * (i - 1) + 1
            w = sw + 2nelemx + 1
            nw = sw + 2 * (2nelemx + 1)

            s = sw + 1
            c = sw + 1 + 2nelemx + 1
            n = sw + 1 + 2 * (2nelemx + 1)

            se = sw + 2
            e = sw + 2 + 2nelemx + 1
            ne = sw + 2 + 2 * (2nelemx + 1)

            element += 1

            element_node[1, element] = sw
            element_node[2, element] = se
            element_node[3, element] = ne
            element_node[4, element] = nw
            element_node[5, element] = s
            element_node[6, element] = e
            element_node[7, element] = n
            element_node[8, element] = w
            element_node[9, element] = c
        end
    end

    return element_node'
end
