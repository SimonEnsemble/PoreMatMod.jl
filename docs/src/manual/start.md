# Getting Started

## Installation

Download and install the [Julia programming language](https://julialang.org/),
v1.6 or higher, and set up [Xtals.jl](https://simonensemble.github.io/Xtals.jl/dev/).

Then, n the Julia REPL, open the package manager (using `]`) and enter:

```julia
pkg> add MOFun
```

Finally, load `MOFun` into the namespace:

```julia
julia> using MOFun # That's it!
```

## Xtals.jl

`MOFun` is built on [`Xtals.jl`](https://github.com/SimonEnsemble/Xtals.jl), which
provides the framework for representing and manipulating crystal structures. It
is recommended that the user be familiar with some of its key functions, which
are re-exported by `MOFun`, particularly [`Crystal`] and [`infer_bonds!`].


## A Guiding Example

This manual makes repeated reference to the following task: given the experimental structure of 
[IRMOF-1](https://raw.githubusercontent.com/SimonEnsemble/MOFun.jl/master/test/data/crystals/IRMOF-1.cif?token=AD3TMGFZCE4WX3J4TDH2BSDAYMO2K), 
produce the hypothetical structure of an isoreticular MOF where the BDC linker has been replaced with a derivative, 2-acetylamido-BDC.
