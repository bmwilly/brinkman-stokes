# A scalable matrix-free Stokes-Brinkman solver in Julia

## Introduction

High-order, parallel, matrix-free finite element solver for Stokes and Brinkman
flow written in Julia.

Details and bechmarks can be found
[here](https://github.com/bmwilly/brinkman-stokes/blob/master/B.%20Williams%20-%20A%20scalable%20matrix-free%20Stokes-Brinkman%20solver%20in%20Julia.pdf).

## Basic usage

There are two implementations of the software, a parallel Julia implementation
and a serial Matlab implementation. For each implementation, there are three
versions that correspond to various finite element assembly methods, using
assembled matrices, matrix-free methods using matrix-vector products, and
matrix-free methods using matrix-matrix products.

Example usage:
```
$ cd julia-parallel/efficient-operators
$ julia
julia> include("sss.jl")
```
will start an interactive session that will prompt
you for a type of test problem to solve (lid-driven cavity or Brinkman flow
through obstacles), mesh size, preconditioning steps, etc.
