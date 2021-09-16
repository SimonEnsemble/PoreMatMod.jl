![logo.JPG](assets/index/logo.JPG)

`PoreMatMod.jl` is a software package in [Julia](https://julialang.org/) for (i) subgraph matching in and (ii) modifying of crystal structures. 
By functioning as a "find-and-replace" tool on atomistic structure models, `PoreMatMod.jl` can search crystals for chemical substructures, create libraries of hypothetical structures, introduce defects into structures, and correct artifacts of X-ray structure determination.

`PoreMatMod.jl` implements 
(i, for find operations) Ullmann's algorithm for subgraph isomorphism search
(ii, for replace operations) the orthogonal Procrustes algorithm for point cloud alignment.  
Periodic boundary conditions are respected, and the unit cell is preserved.

While developed primarily for porous crystals such as metal-organic frameworks (MOFs), `PoreMatMod.jl` can operate on any periodic chemical system as well as discrete molecules.

Example use cases:
:hammer: tuning the chemistry of (e.g., decorating with functional groups) existing crystal structure models to generate libraries of hypothetical materials for computational screening
:hammer: generating heterogeneous, multi-linker MOFs with precise control of functional group placement
:hammer: repairing artifacts in crystal structures---such as missing hydrogen atoms, disorder, and the presence of solvents---determined from X-ray diffraction studies
:hammer: introducing missing-linker and missing-node defects into MOFs to enable computational studies on defect-property relationships
:hammer: searching for subgraphs in libraries of crystal structure models to, e.g., filter structures or characterize the diversity of the library

note: `PoreMatMod.jl` is built on [Xtals.jl](https://github.com/SimonEnsemble/Xtals.jl).

## example: creating a functionalized MOF structure

suppose we wish to decorate the linkers of IRMOF-1 with trifluoromethyl (tfm) groups.
the `PoreMatMod.jl` code below accomplishes this by (i) searching the parent IRMOF-1 structure for a phenylene query fragment and (ii) replacing each instance with a tfm-phenylene replacement fragment to give the child structure.

```julia
# read crystal structure of the parent MOF
parent_xtal = Crystal("IRMOF-1.cif")

# read query and replacement fragments
query_fragment       = moiety("p-phenylene.xyz")  # masked atoms marked with !
replacement_fragment = moiety("tfm-p-phenylene.xyz")

# (1) search parent structure for query fragment
# (2) replace occurrences of query fragment with replacement fragments
#     (with randomly chosen orientations)
child_xtal = replace(parent_xtal, query_fragment => replacement_fragment)
```

![](s_moty-to-r_moty.png)
