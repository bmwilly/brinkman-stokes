# A scalable matrix-free Stokes-Brinkman solver in Julia

## Introduction

High-order, parallel, matrix-free finite element solver for Stokes and Brinkman
flow written in Julia.

Details and benchmarks can be found
[here](https://github.com/bmwilly/brinkman-stokes/blob/master/B.%20Williams%20-%20A%20scalable%20matrix-free%20Stokes-Brinkman%20solver%20in%20Julia.pdf).

## Setup

All commands in this README assume you are in the `julia` REPL in the `BrinkmanStokes` folder:

```shell
cd BrinkmanStokes
julia
```

Install [julia](https://julialang.org/downloads/). Then

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

```julia
julia> ]
pkg> activate .
julia> include("src/sss.jl")
```

will start an interactive session that will prompt
you for a type of test problem to solve (lid-driven cavity or Brinkman flow
through obstacles), mesh size, preconditioning steps, etc.
