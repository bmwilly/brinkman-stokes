"""
Q2-P1 element grid generator for Taylor-Hood finite elements.

Generates pressure node coordinates and element connectivity for Q2-P1
Taylor-Hood elements used in stable Stokes flow discretizations.
"""

using SparseArrays
using Statistics

"""
	q2p1grid(grid::Dict) -> Dict

Generate Q2-P1 element grid from Q2 velocity grid.

Takes a Q2 velocity grid and generates the corresponding pressure grid for
Taylor-Hood Q2-P1 elements. This creates a stable finite element pair for
Stokes flow that satisfies the inf-sup condition.

# Arguments
- `grid::Dict`: Q2 velocity grid containing:
  - `"x"`: x-coordinate vector
  - `"y"`: y-coordinate vector
  - `"xy"`: Velocity node coordinates
  - `"mv"`: Q2 element connectivity
  - `"bound"`: Boundary node indices

# Returns
- `Dict`: Extended grid with pressure information:
  - `"x"`: Updated x-coordinate vector (with corrected mid-points)
  - `"y"`: Updated y-coordinate vector (with corrected mid-points)
  - `"xy"`: Updated velocity node coordinates
  - `"xyp"`: Pressure node coordinates (element centroids)
  - `"ee"`: Element-to-element connectivity (currently empty)

# Q2-P1 Element Structure
- Velocity (Q2): 9 nodes per element (biquadratic)
  - 4 corner nodes + 4 mid-side nodes + 1 center node
- Pressure (P1): 1 node per element (constant per element)
  - Located at element centroid

# Mid-point Correction
For stretched grids, mid-side nodes are repositioned to true edge centers:
- Original: Equally spaced in parameter space
- Corrected: Geometric center of edge endpoints

# Examples
```julia
# Start with Q2 velocity grid
velocity_grid = cavity_domain(3)

# Generate Q2-P1 compatible grid
stokes_grid = q2p1grid(velocity_grid)

# Access pressure nodes
xyp = stokes_grid["xyp"]  # Pressure node coordinates
xy = stokes_grid["xy"]    # Corrected velocity coordinates
```

# Notes
- Pressure nodes are located at element centroids
- Mid-side velocity nodes are corrected for geometric accuracy
- Element connectivity `ee` is currently not implemented (empty array)
- Grid modification ensures geometric consistency

# See also
[`cavity_domain`](@ref), [`stokes_q2p1`](@ref)
"""
function q2p1grid(grid::Dict)
    # Extract velocity grid data
    x = copy(grid["x"])
    y = copy(grid["y"])
    xy = copy(grid["xy"])
    mv = grid["mv"]
    bound = grid["bound"]

    # Grid dimensions
    nvtx = size(xy, 1)        # Number of velocity nodes
    nel = size(mv, 1)         # Number of elements
    nx, ny = length(x), length(y)

    println("Generating Q2-P1 grid from Q2 velocity grid")
    println("  Velocity nodes: $(nvtx)")
    println("  Elements: $(nel)")
    println("  Grid size: $(nx)×$(ny)")

    # Correct mid-side points for stretched grids
    println("  Correcting mid-side node positions...")

    # Extract coordinate arrays for easier manipulation
    xx = xy[:, 1]
    yy = xy[:, 2]

    # Correct y-direction mid-side points
    for k in 2:2:ny
        y_old = y[k]
        y_new = 0.5 * (y[k + 1] + y[k - 1])  # True geometric center

        # Find and update all nodes with this y-coordinate
        node_indices = findall(yy .≈ y_old)
        yy[node_indices] .= y_new
        y[k] = y_new
    end

    # Correct x-direction mid-side points
    for k in 2:2:nx
        x_old = x[k]
        x_new = 0.5 * (x[k + 1] + x[k - 1])  # True geometric center

        # Find and update all nodes with this x-coordinate
        node_indices = findall(xx .≈ x_old)
        xx[node_indices] .= x_new
        x[k] = x_new
    end

    # Update coordinate matrix with corrected positions
    xy = [xx yy]

    # Generate pressure node coordinates (element centroids)
    println("  Computing pressure node coordinates...")

    xc = zeros(nel)
    yc = zeros(nel)

    for ielem in 1:nel
        # Get corner node coordinates for this element
        corner_nodes = mv[ielem, 1:4]  # Q2 corner nodes

        # Compute centroid of corner nodes
        xc[ielem] = mean(xx[corner_nodes])
        yc[ielem] = mean(yy[corner_nodes])
    end

    xyp = [xc yc]  # Pressure node coordinate matrix

    # Element-to-element connectivity (currently not implemented)
    # This would be useful for advanced algorithms but is computationally expensive
    # to compute and not currently needed

    println("  Computing element connectivity...")
    ee = Int[]  # Placeholder - could implement adjacency matrix if needed

    # The following commented code shows how element connectivity could be computed:
    # adj = spzeros(nvtx, nvtx)
    #
    # # Build adjacency matrix based on shared edges
    # for el in 1:nel
    #     nodes = mv[el, :]
    #     # Add connections for each edge
    #     adj[nodes[1], nodes[2]] = el  # Bottom edge
    #     adj[nodes[2], nodes[3]] = el  # Right edge
    #     adj[nodes[3], nodes[4]] = el  # Top edge
    #     adj[nodes[4], nodes[1]] = el  # Left edge
    # end
    #
    # # Extract element neighbors from adjacency matrix
    # ee = zeros(Int, nel, 4)  # 4 neighbors per element
    # for el in 1:nel
    #     # Find neighboring elements for each edge
    #     # Implementation would go here...
    # end

    println("Q2-P1 grid generation complete")
    println("  Pressure nodes: $(nel)")
    println("  Mid-side corrections applied")

    # Return updated grid with pressure information
    return Dict(
        "x" => x,      # Corrected x-coordinate vector
        "y" => y,      # Corrected y-coordinate vector
        "xy" => xy,    # Corrected velocity node coordinates
        "xyp" => xyp,  # Pressure node coordinates
        "ee" => ee,     # Element connectivity (empty)
    )
end
