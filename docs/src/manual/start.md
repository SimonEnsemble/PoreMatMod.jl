# Getting Started

## Installation and Use

Download and install the [Julia programming language](https://julialang.org/),
v1.5 or higher.

In the Julia REPL, open the package manager (using `]`) and enter:

```julia
pkg> add MOFun
```

Next, load `MOFun` into the namespace:

```julia
julia> using MOFun # That's it!
```

## Xtals.jl

`MOFun` is built on [`Xtals.jl`](https://github.com/SimonEnsemble/Xtals.jl), which
provides the framework for representing and manipulating crystal structures. It
is recommended that the user be familiar with some of its key functions, which
are re-exported by `MOFun`:

- [`Xtals.set_path_to_data`](@ref)

- [`Xtals.Crystal`](@ref)

- [`Xtals.infer_bonds!`](@ref)

- [`Xtals.infer_geometry_based_bonds!`](@ref)

- [`Xtals.write_cif`](@ref)

- [`Xtals.write_vtk`](@ref)

- [`Xtals.write_xyz`](@ref)

- [`Xtals.write_bond_information`](@ref)

## A Guiding Example

This manual makes repeated reference to the following task: given the experimental
structure of [IRMOF-1](assets/IRMOF-1.cif), produce the hypothetical structure of
an isoreticular MOF where the BDC linker has been replaced with a derivative,
2-acetylamido-BDC.
