# MOFun.jl

A pure-[Julia](https://julialang.org/) software package for manipulating chemical
structures based on an open-source substructure find-and-replace algorithm.

![logo]()

`MOFun.jl` can identify chemical substructures, create hypothetical structure
libraries, and repair common defects in experimental data.  For example, it can
clean and derivatize

![example case: disordered w/ guest -> activated derivative]()

```julia
using MOFun
xtal = Crystal("guest&disorder.cif", infer_bonds=:periodic)
repaired = (moiety("disordered!.xyz") => moiety("discrete.xyz")) ∈ xtal
active = (moiety("solvent.xyz") => moiety("nothing.xyz")) ∈ repaired
novel = (moiety("2-!-p-phenylene") => moiety("2-Me-p-phenylene")) ∈ active
```
