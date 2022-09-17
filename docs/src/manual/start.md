# Getting Started

## Installation

Download and install the [Julia programming language](https://julialang.org/), at version 1.6 or higher.
`PoreMatMod.jl` is not currently stable on version 1.7.

To install `PoreMadMod.jl` (officially registered as a Julia package), in the Julia REPL, enter the package manager by typing `]` and enter:

```
pkg> add PoreMatMod
```

## Loading the `PoreMatMod.jl` package

To load the `PoreMatMod.jl` package, so that its functions are imported into your namespace, in your Julia code:

```julia
julia> using PoreMatMod
```

We recommend writing Julia code and performing find-and-replace tasks with `PoreMadMod.jl` using interactive [Pluto notebooks](https://github.com/fonsp/Pluto.jl).

## Running tests (optional)

Run the unit tests associated with `PoreMatMod.jl` by entering package mode in the Julia REPL via (`]`) and entering:

```
pkg> test PoreMatMod
```

!!! note

    `PoreMatMod.jl` is built on [`Xtals.jl`](https://github.com/SimonEnsemble/Xtals.jl), which provides:
    * the data structure and reader, `Crystal`, for crystal structures 
    * the `infer_bonds!` function that assigns bonds between atoms of a `Crystal`
