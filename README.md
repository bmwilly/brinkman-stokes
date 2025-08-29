# A scalable matrix-free Stokes-Brinkman solver in Julia

## Introduction

High-order, parallel, matrix-free finite element solver for Stokes and Brinkman
flow written in Julia.

Details and benchmarks can be found
[here](https://github.com/bmwilly/brinkman-stokes/blob/master/B.%20Williams%20-%20A%20scalable%20matrix-free%20Stokes-Brinkman%20solver%20in%20Julia.pdf).

**NOTE**: This code and paper was for my master's thesis. It has been updated to use modern Julia. To view the code as it existed when I submitted my thesis, see the [thesis-submission](https://github.com/bmwilly/brinkman-stokes/tree/thesis-submission) tag.

## Setup

### Install julia

Install julia using the [official instructions](<https://julialang.org/downloads/>), which will also install `juliaup`. Then

```shell
juliaup add 1.11
juliaup default 1.11
```

### Install Python

Python is required for the plotting scripts. Install Python using your favorite method. We recommend using [uv](https://docs.astral.sh/uv/):

```shell
curl -LsSf https://astral.sh/uv/install.sh | sh  # install uv
uv python install 3.12
uv venv .venv
source .venv/bin/activate
uv sync
```

Then in julia:

```julia
ENV["PYTHON"] = "/Users/bwilliams/projects/personal/brinkman-stokes/.venv/bin/python"  # or wherever your python is

using Pkg

Pkg.add("PyCall")
Pkg.build("PyCall")
```

### Install packages

```shell
julia --project=BrinkmanStokes
```

Then

```julia
julia> ]
(BrinkmanStokes) pkg> instantiate
(BrinkmanStokes) pkg> precompile
```

#### Updating packages

Use `juliaup` to install and use the desired version of Julia:

```shell
juliaup add <new_version>
juliaup default <new_version>
```

Then start Julia:

```shell
julia --project=BrinkmanStokes
```

and update the packages:

```julia
julia> ]
(BrinkmanStokes) pkg> update
(BrinkmanStokes) pkg> instantiate
(BrinkmanStokes) pkg> resolve
```

Then commit the changes to the `Manifest.toml` and `Project.toml` files.

## Code Formatting

This project uses [Runic.jl](https://github.com/fredrikekre/Runic.jl) for consistent code formatting across all Julia files.

### Installing Runic.jl

Since Runic.jl is not available in the General registry, it needs to be installed in a separate environment:

```shell
julia --project=@runic -e 'using Pkg; Pkg.add(url = "https://github.com/fredrikekre/Runic.jl")'
```

### Formatting Code

To format all Julia files in the project:

```shell
julia --project=@runic -e 'using Runic; exit(Runic.main(["--inplace", "."]))'
```

This will recursively format all `.jl` files in the project directory using Runic's fixed formatting rules.

### Optional: Shell Alias

For convenience, you can add this alias to your shell configuration (`.bashrc`, `.zshrc`, etc.):

```shell
alias runic="julia --project=@runic -e 'using Runic; exit(Runic.main(ARGS))'"
```

Then you can simply run:

```shell
runic --inplace .
```

## Basic usage

There are two implementations of the software, a parallel Julia implementation
and a serial Matlab implementation. For each implementation, there are three
versions that correspond to various finite element assembly methods, using
assembled matrices, matrix-free methods using matrix-vector products, and
matrix-free methods using matrix-matrix products.

Example usage:

```shell
julia --project=BrinkmanStokes
```

```julia
julia> include("BrinkmanStokes/src/sss.jl")
```

will start an interactive session that will prompt
you for a type of test problem to solve (lid-driven cavity or Brinkman flow
through obstacles), mesh size, preconditioning steps, etc.
