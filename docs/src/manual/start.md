# Getting Started

## Installation

Download and install the [Julia programming language](https://julialang.org/), v1.6 or higher.
Then, in the Julia REPL, open the package manager (using `]`) and enter:

```julia
pkg> add PoreMatMod
```

Finally, load `PoreMatMod.jl` into the namespace:

```julia
julia> using PoreMatMod # That's it!
```

## Xtals.jl

`PoreMatMod.jl` is built on [`Xtals.jl`](https://github.com/SimonEnsemble/Xtals.jl), which provides the framework for representing and manipulating crystal structures. 
It is recommended that the user be familiar with some of its key functions, which are re-exported by `PoreMatMod`, particularly [`Crystal`] and [`infer_bonds!`].


## A Guiding Example

This manual makes repeated reference to the following task: given the experimentally-determined structure of [IRMOF-1](../../../assets/start/IRMOF-1.cif), produce the hypothetical structure of an isoreticular MOF where the BDC linker has been replaced with a derivative, 2-acetylamido-BDC.

![IRMOF-1_example](../../assets/start/example1.png)
