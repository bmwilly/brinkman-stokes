# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Repository Overview

This is a scalable matrix-free Stokes-Brinkman solver implemented in both Julia and MATLAB. The code is a modern Julia adaptation of a master's thesis project that implements high-order, parallel, matrix-free finite element methods for Stokes and Brinkman flow problems.

## Quick Start

### Julia Development (Primary Implementation)

**Setup Environment:**

```bash
julia --project=BrinkmanStokes
```

**Install Dependencies:**

```julia
julia> ]
(BrinkmanStokes) pkg> instantiate
(BrinkmanStokes) pkg> precompile
```

**Run Interactive Solver:**

```julia
julia> include("BrinkmanStokes/src/sss.jl")
```

**Run Tests:**

```julia
julia> include("BrinkmanStokes/src/test.jl")
```

### Python Environment (For Plotting)

**Setup Python Environment:**

```bash
uv venv .venv
source .venv/bin/activate
uv sync
```

**Configure PyCall in Julia:**

```julia
ENV["PYTHON"] = "/path/to/brinkman-stokes/.venv/bin/python"
using Pkg; Pkg.add("PyCall"); Pkg.build("PyCall")
```

## Architecture

### Core Structure

- **BrinkmanStokes/**: Main Julia package containing the solver implementation
- **matlab/**: Three MATLAB implementations for comparison (global-matrix, matrix-free-stacked, matrix-free-tensor)
- **BrinkmanStokes/src/homg/**: High-order finite element framework (Julia port of HoMG)
- **BrinkmanStokes/src/efficient-operators/**: Core solver algorithms

### Key Components

1. **Solver Entry Point** (`sss.jl`): Interactive session that prompts for problem parameters
2. **Domain Types**:
   - Domain 1: Lid-driven cavity flow (Stokes)
   - Domain 2: Brinkman flow through obstacles
3. **Matrix-Free Methods**: Three implementation approaches for finite element assembly
4. **Geometric Multigrid**: Preconditioning using `mg_diff.jl` and `m_st_mg.jl`
5. **Output System**: Unified file organization via `output_utils.jl`

### Finite Element Framework (homg)

The `homg/` directory contains a high-order finite element framework supporting:

- Hexahedral elements with arbitrary polynomial degree
- Tensor product basis functions
- Coordinate transformations
- Grid generation utilities

### Solver Algorithms

- **GMRES** iterative solver with optional geometric multigrid preconditioning
- Matrix-free operators for memory efficiency
- Convergence monitoring and plotting
- Support for both assembled and matrix-free approaches

## Development Commands

### Package Management

**Update Julia packages:**

```bash
julia --project=BrinkmanStokes -e 'using Pkg; Pkg.update(); Pkg.instantiate(); Pkg.resolve()'
```

### Running Solvers

**Interactive solver with user input:**

```julia
include("BrinkmanStokes/src/sss.jl")
```

**Direct solve (programmatic):**

```julia
using BrinkmanStokes
include("BrinkmanStokes/src/efficient-operators/stokes_flow/solve_stokes.jl")
sol = solve_stokes(1, 16)  # domain=1 (cavity), msize=16
```

### Output Management

All solver outputs use a standardized directory structure:

- `BrinkmanStokes/output/{implementation}/domain={domain}/size={msize}/`
- Convergence plots saved as PNG
- Convergence data saved as CSV
- Use `get_output_file()` and `get_output_dir()` functions for consistency

## Key Files

- `BrinkmanStokes/src/sss.jl`: Main interactive entry point
- `BrinkmanStokes/src/efficient-operators/stokes_flow/solve_stokes.jl`: Core solver
- `BrinkmanStokes/src/output_utils.jl`: Output file management utilities
- `BrinkmanStokes/Project.toml`: Julia package dependencies
- `pyproject.toml`: Python dependencies for plotting

## MATLAB Implementations

Three parallel MATLAB implementations exist for comparison:

1. **global-matrix/**: Traditional assembled matrix approach
2. **matrix-free-stacked/**: Matrix-free using stacked operations
3. **matrix-free-tensor/**: Matrix-free using tensor products

Each MATLAB version has its own `solve_square_stokes.m` entry point.
