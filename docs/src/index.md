# MOFun.jl

![logo]()

A pure-[Julia](https://julialang.org/) software package for manipulating chemical
structures of porous crystals.  Built on the Simon Ensemble
[Xtals](https://github.com/SimonEnsemble/Xtals.jl) software to clean experimental
and calculated data and provide novel hypothetical structural inputs to
[PorousMaterials](https://github.com/SimonEnsemble/PorousMaterials.jl).  It is
intended primarily for MOFs and other porous crystalline materials, but works
with discrete molecular structures and simulated ensembles as well.

`MOFun.jl` can identify chemical substructures, create hypothetical structure
libraries, and repair common defects in experimental data, using an implementation
of Ullmann's algorithm for substructure isomorphism and the orthogonal Procrustes
algorithm for point cloud alignment.

Example: repairing, activating, and functionalizing an experimental structure.

![example case: disordered w/ guest -> activated derivative]()
(link to pub for CCDC EMEHUB)

```julia
using MOFun
xtal = Crystal("EMEHUB.cif", infer_bonds=true)
repaired = (moiety("disordered!") => moiety("discrete")) ∈ xtal
active = (moiety("solvent") => moiety(nothing)) ∈ repaired
novel = (moiety("2-!-p-phenylene") => moiety("2-Me-p-phenylene")) ∈ active
```
