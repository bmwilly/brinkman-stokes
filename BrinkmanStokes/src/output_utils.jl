"""
Output utilities for consistent file organization across implementations
"""

"""
	get_output_dir(implementation, domain, msize; subdir="")

Generate standardized output directory path for the given implementation.

# Arguments
- `implementation::String`: One of "efficient-operators", "matrix-free-stacked", "matrix-free-tensor"
- `domain::Int`: Domain type (1=lid-driven cavity, 2=brinkman)
- `msize::Int`: Mesh size parameter
- `subdir::String`: Optional subdirectory (e.g., "plots")

# Returns
- `String`: Full path to output directory
"""
function get_output_dir(implementation::String, domain::Int, msize::Int; subdir::String = "")
	# Get BrinkmanStokes root directory (assumes this file is in src/)
	repo_root = dirname(@__DIR__)

	# Build path: BrinkmanStokes/output/{implementation}/domain={domain}/size={msize}/{subdir}
	output_path = joinpath(repo_root, "output", implementation, "domain=$domain", "size=$msize")

	if !isempty(subdir)
		output_path = joinpath(output_path, subdir)
	end

	# Ensure directory exists
	mkpath(output_path)

	return output_path
end

"""
	get_output_file(implementation, domain, msize, filename; subdir="")

Generate standardized output file path for the given implementation.

# Arguments
- `implementation::String`: One of "efficient-operators", "matrix-free-stacked", "matrix-free-tensor"
- `domain::Int`: Domain type (1=lid-driven cavity, 2=brinkman)
- `msize::Int`: Mesh size parameter
- `filename::String`: Name of the output file
- `subdir::String`: Optional subdirectory (e.g., "plots")

# Returns
- `String`: Full path to output file
"""
function get_output_file(implementation::String, domain::Int, msize::Int, filename::String; subdir::String = "")
	output_dir = get_output_dir(implementation, domain, msize; subdir = subdir)
	return joinpath(output_dir, filename)
end
