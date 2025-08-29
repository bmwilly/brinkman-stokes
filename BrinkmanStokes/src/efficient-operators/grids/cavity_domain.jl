"""
Square cavity Q2 grid generator.

Generates structured quadrilateral meshes for square cavity domains using
biquadratic (Q2) elements. The domain is [-1,1] × [-1,1] with symmetric
grid refinement.
"""

include("../helpers/input.jl")
include("../helpers/meshgrid.jl")

"""
	cavity_domain(msize::Int) -> Dict

Generate Q2 finite element grid for square cavity domain.

Creates a structured mesh of the unit square domain [-1,1] × [-1,1] using
biquadratic (Q2) elements with 9 nodes per element. The mesh has symmetric
refinement and power-of-2 structure suitable for multigrid methods.

# Arguments
- `msize::Int`: Mesh size parameter (≥ 1). Grid has 2^(msize+1) intervals per direction

# Returns
- `Dict`: Dictionary containing complete grid information:
  - `"mv"`: Element connectivity matrix (nel × 9)
  - `"xy"`: Node coordinates (nnodes × 2)
  - `"bound"`: Boundary node indices
  - `"mbound"`: Boundary element information
  - `"x"`, `"y"`: Coordinate vectors for each direction
  - `"bndxy"`: Boundary vertex coordinates for domain outline
  - `"bnde"`: Boundary edge connectivity for graphics
  - `"obs"`: Obstacle information (empty for square cavity)

# Grid Structure
- Domain: [-1,1] × [-1,1]
- Elements: np² where np = 2^msize
- Nodes per element: 9 (Q2 biquadratic)
- Node ordering: corners (1-4), mid-sides (5-8), center (9)

# Boundary Identification
Boundaries are labeled as:
1. Bottom edge: y = -1
2. Right edge: x = +1
3. Top edge: y = +1
4. Left edge: x = -1

# Examples
```julia
# Generate coarse 4×4 element mesh
grid = cavity_domain(2)  # 2^2 = 4 elements per direction

# Generate fine 16×16 element mesh
grid = cavity_domain(4)  # 2^4 = 16 elements per direction

# Access grid data
mv = grid["mv"]        # Element connectivity
xy = grid["xy"]        # Node coordinates
bound = grid["bound"]  # Boundary nodes
```

# Notes
- Grid is symmetric about origin with uniform spacing
- Compatible with geometric multigrid coarsening
- Mid-side nodes are placed at element edge centers
- Suitable for Taylor-Hood Q2-P1 Stokes elements

# See also
[`q2p1grid`](@ref), [`meshgrid`](@ref)
"""
function cavity_domain(msize::Int)
	# Input validation
	msize >= 1 || throw(ArgumentError("msize must be ≥ 1"))
	msize <= 8 || @warn "Large msize ($(msize)) may require significant memory"

	println("Generating Q2 grid for square cavity domain")

	# Calculate grid parameters
	n = 2^(msize + 1)           # Total intervals per direction
	np = Int(n / 2)             # Elements per direction
	nel = np^2                  # Total number of elements

	println("  Elements per direction: $(np)")
	println("  Total elements: $(nel)")
	println("  Grid intervals: $(n) per direction")

	# Generate coordinate vectors
	# Symmetric grid: [-1, -1+h, ..., -h, 0, h, ..., 1-h, 1]
	h = 1.0 / np                # Element size
	yy = collect(h:h:1.0)       # Positive y coordinates
	ypos = [0.0; yy]            # Include origin
	yneg = -reverse(yy)         # Negative y coordinates
	y = [yneg; ypos]            # Complete y vector
	x = copy(y)                 # Square domain: x = y

	# Generate nodal coordinates using meshgrid
	nvtx = (n + 1)^2           # Total number of nodes
	X, Y = meshgrid(x, y)

	# Flatten coordinate matrices to vectors
	xx = reshape(X', nvtx)
	yy = reshape(Y', nvtx)
	xy = [xx yy]               # Nodal coordinate matrix

	# Generate element connectivity matrix
	println("  Generating element connectivity...")

	mv = zeros(Int64, nel, 9)  # 9 nodes per Q2 element
	kx = 1                     # Current x index
	ky = 1                     # Current y index
	mel = 0                    # Element counter

	# Loop over elements in row-major order
	for j in 1:np              # y-direction (rows)
		for i in 1:np          # x-direction (columns)
			# Reference node (bottom-left corner of element)
			mref = (n + 1) * (ky - 1) + kx
			mel += 1

			# Q2 element node numbering (local):
			#  4 --- 7 --- 3
			#  |     |     |
			#  8 --- 9 --- 6
			#  |     |     |
			#  1 --- 5 --- 2

			# Map local to global node numbers
			nvv = zeros(Int, 9)
			nvv[1] = mref                    # Bottom-left corner
			nvv[2] = mref + 2                # Bottom-right corner
			nvv[3] = mref + 2n + 4           # Top-right corner
			nvv[4] = mref + 2n + 2           # Top-left corner
			nvv[5] = mref + 1                # Bottom edge center
			nvv[6] = mref + n + 3            # Right edge center
			nvv[7] = mref + 2n + 3           # Top edge center
			nvv[8] = mref + n + 1            # Left edge center
			nvv[9] = mref + n + 2            # Element center

			mv[mel, :] = nvv
			kx += 2                          # Move to next element in x
		end
		ky += 2                              # Move to next row
		kx = 1                               # Reset x position
	end

	# Identify boundary nodes and elements
	println("  Identifying boundary nodes...")

	# Bottom boundary (y = -1)
	k1 = findall(xy[:, 2] .== -1.0)
	e1 = Int[]
	for k in 1:nel
		if any(in(k1), mv[k, 5])  # Check if mid-side node 5 is on boundary
			push!(e1, k)
		end
	end
	ef1 = ones(Int, length(e1))  # Boundary flag = 1

	# Right boundary (x = +1, excluding corners)
	k2 = findall((xy[:, 1] .== 1.0) .& (-1.0 .< xy[:, 2]) .& (xy[:, 2] .< 1.0))
	e2 = Int[]
	for k in 1:nel
		if any(in(k2), mv[k, 6])  # Check if mid-side node 6 is on boundary
			push!(e2, k)
		end
	end
	ef2 = 2 * ones(Int, length(e2))  # Boundary flag = 2

	# Top boundary (y = +1)
	k3 = findall(xy[:, 2] .== 1.0)
	e3 = Int[]
	for k in 1:nel
		if any(in(k3), mv[k, 7])  # Check if mid-side node 7 is on boundary
			push!(e3, k)
		end
	end
	ef3 = 3 * ones(Int, length(e3))  # Boundary flag = 3

	# Left boundary (x = -1, excluding corners)
	k4 = findall((xy[:, 1] .== -1.0) .& (-1.0 .< xy[:, 2]) .& (xy[:, 2] .< 1.0))
	e4 = Int[]
	for k in 1:nel
		if any(in(k4), mv[k, 8])  # Check if mid-side node 8 is on boundary
			push!(e4, k)
		end
	end
	ef4 = 4 * ones(Int, length(e4))  # Boundary flag = 4

	# Collect all boundary nodes
	bound = sort([k1; k2; k3; k4])
	mbound = [e1' ef1'; e2' ef2'; e3' ef3'; e4' ef4']

	# Define domain boundary for visualization
	# Square domain vertices in counterclockwise order
	bndxy = [
		-1.0 -1.0;   # Bottom-left
		1.0 -1.0;   # Bottom-right
		1.0  1.0;   # Top-right
		-1.0  1.0
	]   # Top-left

	# Boundary edges for graphics (node1, node2, boundary_id)
	bnde = [
		1 2 1;        # Bottom edge
		2 3 2;        # Right edge
		3 4 3;        # Top edge
		4 1 4
	]        # Left edge

	# No obstacles in square cavity
	obs = Int[]

	# All edges may need stretching consideration
	sbnde = [1, 2, 3, 4]

	println("Grid generation complete:")
	println("  Total nodes: $(size(xy, 1))")
	println("  Boundary nodes: $(length(bound))")
	println("  Elements: $(nel)")

	# Return complete grid data structure
	return Dict(
		"mv" => mv,           # Element connectivity matrix
		"xy" => xy,           # Node coordinates
		"bound" => bound,     # Boundary node indices
		"mbound" => mbound,   # Boundary element information
		"x" => x,             # x-coordinate vector
		"y" => y,             # y-coordinate vector
		"bndxy" => bndxy,     # Domain boundary vertices
		"bnde" => bnde,       # Boundary edge connectivity
		"obs" => obs,          # Obstacle information (empty)
	)
end
