# MOFun.jl

A pure-[Julia](https://julialang.org/) software package for manipulating chemical
structures of porous crystals.  `MOFun` is built on the
[@SimonEnsemble](https://SimonEnsemble.github.io) software
[Xtals](https://github.com/SimonEnsemble/Xtals.jl) to clean experimental and
calculated data and provide novel hypothetical structural inputs to
[PorousMaterials](https://github.com/SimonEnsemble/PorousMaterials.jl).  It is
intended primarily for MOFs and other porous crystalline materials, but can work
with other periodic structures, ensembles, and discrete molecules as well.

`MOFun.jl` can identify chemical substructures, create hypothetical structure
libraries, and correct disorder in experimental data, using an implementation
of Ullmann's algorithm for substructure isomorphism and the orthogonal Procrustes
algorithm for point cloud alignment.  Periodic cell boundaries are treated
automatically, and the unit cell is preserved in transformations.

**Example**: repairing, activating, and functionalizing an experimental
structure.  [This structure](https://dx.doi.org/10.5517/ccdc.csd.cc1ldj8s), of
a MOF called [SIFSIX-1-Cu](https://dx.doi.org/10.1126/science.aaf2458), contains
disordered PyC2 linkers and acetylene guest molecules.

![example]()

Loading the data, resolving the disorder, removing the guest molecules, replacing
the linkers, and saving the result can be done with a very short script:

```julia
# Import the module
using MOFun
# Load some messy data
xtal = Crystal("EMEHUB.cif", infer_bonds=:voronoi, periodic_bonds=true)
# Repair the disordered linkers
repaired = (moiety("disordered_bipy!") => moiety("discrete")) ∈ xtal
# Remove the guest molecules
active = (moiety("acetylene") => moiety(nothing)) ∈ repaired
# Add a functional group
novel = (moiety("2-H!-PyC2") => moiety("2-Me-PyC2")) ∈ active
# Save the result
write_cif(novel, "2-Me-SIFSIX-1-Cu")
```
