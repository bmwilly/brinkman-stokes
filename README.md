# A scalable matrix-free Stokes-Brinkman solver in Julia

## Introduction

High-order, parallel, matrix-free finite element solver for Stokes and Brinkman
flow written in Julia.

Details and benchmarks can be found
[here](https://github.com/bmwilly/brinkman-stokes/blob/master/B.%20Williams%20-%20A%20scalable%20matrix-free%20Stokes-Brinkman%20solver%20in%20Julia.pdf).

## Setup

### Install julia

Install julia using the [official instructions](<https://julialang.org/downloads/>), which will also install `juliaup`. Then

```shell
juliaup add 1.10.2
juliaup default 1.10.2
```

### Install Python

Python is required for the plotting scripts. Install Python using your favorite method. We recommend using `pyenv`:

```shell
brew install pyenv  # on macOS
pyenv install 3.12.3
pyenv local 3.12.3
pyenv -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Then in julia:

```julia
ENV["PYTHON"] = "/Users/bwilliams/projects/personal/brinkman-stokes/.venv/bin/python"  # or wherever your python is

using Pkg

Pkg.add("PyCall")
Pkg.build("PyCall")
```

### Install package and dependencies

```shell
cd BrinkmanStokes
julia
```

Then

```julia
julia> ]
pkg> activate .
(BrinkmanStokes) pkg> instantiate
(BrinkmanStokes) pkg> precompile
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
