```@meta
DocTestSetup = quote
    using PoreMatMod
    parent = Crystal("IRMOF-1.cif")
    infer_bonds!(parent, true)
    query = moiety("p-phenylene.xyz")
end
```

# Subgraph matching (substructure searches)

![find graphic](../../assets/find/s_moty-in-xtal.png)

`PoreMatMod.jl` conducts subgraph matching, i.e. searches for subgraphs of a `parent` graph isomorphic to a `query` graph, using [Ullmann's algorithm for subgraph isomorphisms](https://doi.org/10.1145/321921.321925).

During subgraph matching, both the `parent` crystal structure and `query` fragment are represented by node-labeled (by the chemical species) graphs. For crystals, bonds across the unit cell boundaries of periodic materials are accounted for, allowing us to find subgraph isomorphisms when the fragment is split across a unit cell boundary.

# Substructure Searches: how to

To learn by example, suppose we wish to search the IRMOF-1 crystal structure for *p*-phenylene fragments.

First, we load the `query` fragment and `parent` structure:
```julia
parent = Crystal("IRMOF-1.cif")
infer_bonds!(parent_xtal, true) # true to infer bonds across the periodic boundary

query = moiety("p-phenylene.xyz")
```

With a `parent` and `query` loaded, execute a search:

```jldoctest find
search = substructure_search(query, parent)
# output
p-phenylene.xyz ∈ IRMOF-1.cif
96 hits in 24 locations.
```

!!! note "Syntactic sugar for substructure search"
    The `∈` (`\in` then hit `Tab` for this Unicode character) infix operator will also execute the search:

    ```julia
    search = query ∈ parent
    ```


Both functions `substructure_search` and `∈`` return a `Search` object which stores the `query` and `parent` structures and the result of the search, `isomorphisms`, a nested vector giving the query-to-parent correpondences.

```jldoctest find
search.isomorphisms
# output
24-element Vector{Vector{Vector{Int64}}}:
 [[233, 306, 318, 245, 185, 197, 414, 329, 402, 341], [245, 318, 306, 233, 197, 185, 402, 341, 414, 329], [306, 233, 245, 318, 185, 197, 341, 402, 329, 414], [318, 245, 233, 306, 197, 185, 329, 414, 341, 402]]
 [[234, 305, 317, 246, 186, 198, 413, 330, 401, 342], [246, 317, 305, 234, 198, 186, 401, 342, 413, 330], [305, 234, 246, 317, 186, 198, 342, 401, 330, 413], [317, 246, 234, 305, 198, 186, 330, 413, 342, 401]]
 [[235, 308, 319, 248, 187, 200, 415, 331, 404, 344], [248, 319, 308, 235, 200, 187, 404, 344, 415, 331], [308, 235, 248, 319, 187, 200, 344, 404, 331, 415], [319, 248, 235, 308, 200, 187, 331, 415, 344, 404]]
 [[236, 307, 320, 247, 188, 199, 416, 332, 403, 343], [247, 320, 307, 236, 199, 188, 403, 343, 416, 332], [307, 236, 247, 320, 188, 199, 343, 403, 332, 416], [320, 247, 236, 307, 199, 188, 332, 416, 343, 403]]
 [[237, 262, 280, 255, 189, 207, 376, 333, 358, 351], [255, 280, 262, 237, 207, 189, 358, 351, 376, 333], [262, 237, 255, 280, 189, 207, 351, 358, 333, 376], [280, 255, 237, 262, 207, 189, 333, 376, 351, 358]]
 [[238, 261, 279, 256, 190, 208, 375, 334, 357, 352], [256, 279, 261, 238, 208, 190, 357, 352, 375, 334], [261, 238, 256, 279, 190, 208, 352, 357, 334, 375], [279, 256, 238, 261, 208, 190, 334, 375, 352, 357]]
 [[239, 264, 277, 254, 191, 206, 373, 335, 360, 350], [254, 277, 264, 239, 206, 191, 360, 350, 373, 335], [264, 239, 254, 277, 191, 206, 350, 360, 335, 373], [277, 254, 239, 264, 206, 191, 335, 373, 350, 360]]
 [[240, 263, 278, 253, 192, 205, 374, 336, 359, 349], [253, 278, 263, 240, 205, 192, 359, 349, 374, 336], [263, 240, 253, 278, 192, 205, 349, 359, 336, 374], [278, 253, 240, 263, 205, 192, 336, 374, 349, 359]]
 [[241, 290, 299, 252, 193, 204, 395, 337, 386, 348], [252, 299, 290, 241, 204, 193, 386, 348, 395, 337], [290, 241, 252, 299, 193, 204, 348, 386, 337, 395], [299, 252, 241, 290, 204, 193, 337, 395, 348, 386]]
 [[242, 289, 300, 251, 194, 203, 396, 338, 385, 347], [251, 300, 289, 242, 203, 194, 385, 347, 396, 338], [289, 242, 251, 300, 194, 203, 347, 385, 338, 396], [300, 251, 242, 289, 203, 194, 338, 396, 347, 385]]
 ⋮
 [[260, 283, 296, 271, 212, 219, 392, 356, 379, 367], [271, 296, 283, 260, 219, 212, 379, 367, 392, 356], [283, 260, 271, 296, 212, 219, 367, 379, 356, 392], [296, 271, 260, 283, 219, 212, 356, 392, 367, 379]]
 [[265, 314, 323, 276, 213, 224, 419, 361, 410, 372], [276, 323, 314, 265, 224, 213, 410, 372, 419, 361], [314, 265, 276, 323, 213, 224, 372, 410, 361, 419], [323, 276, 265, 314, 224, 213, 361, 419, 372, 410]]
 [[266, 313, 324, 275, 214, 223, 420, 362, 409, 371], [275, 324, 313, 266, 223, 214, 409, 371, 420, 362], [313, 266, 275, 324, 214, 223, 371, 409, 362, 420], [324, 275, 266, 313, 223, 214, 362, 420, 371, 409]]
 [[267, 316, 322, 273, 215, 221, 418, 363, 412, 369], [273, 322, 316, 267, 221, 215, 412, 369, 418, 363], [316, 267, 273, 322, 215, 221, 369, 412, 363, 418], [322, 273, 267, 316, 221, 215, 363, 418, 369, 412]]
 [[268, 315, 321, 274, 216, 222, 417, 364, 411, 370], [274, 321, 315, 268, 222, 216, 411, 370, 417, 364], [315, 268, 274, 321, 216, 222, 370, 411, 364, 417], [321, 274, 268, 315, 222, 216, 364, 417, 370, 411]]
 [[285, 310, 328, 303, 225, 231, 424, 381, 406, 399], [303, 328, 310, 285, 231, 225, 406, 399, 424, 381], [310, 285, 303, 328, 225, 231, 399, 406, 381, 424], [328, 303, 285, 310, 231, 225, 381, 424, 399, 406]]
 [[286, 309, 327, 304, 226, 232, 423, 382, 405, 400], [304, 327, 309, 286, 232, 226, 405, 400, 423, 382], [309, 286, 304, 327, 226, 232, 400, 405, 382, 423], [327, 304, 286, 309, 232, 226, 382, 423, 400, 405]]
 [[287, 312, 325, 302, 227, 230, 421, 383, 408, 398], [302, 325, 312, 287, 230, 227, 408, 398, 421, 383], [312, 287, 302, 325, 227, 230, 398, 408, 383, 421], [325, 302, 287, 312, 230, 227, 383, 421, 398, 408]]
 [[288, 311, 326, 301, 228, 229, 422, 384, 407, 397], [301, 326, 311, 288, 229, 228, 407, 397, 422, 384], [311, 288, 301, 326, 228, 229, 397, 407, 384, 422], [326, 301, 288, 311, 229, 228, 384, 422, 397, 407]]
```

In this example, the `query` fragment (*p*-phenylene) occurs in 24 different location in the `parent` crystal structure, with 4 symmetry-equivalent isomorphisms at each location, for a total of 96 subgraph isomorphisms.

The number of locations---the number of unique substructures of the `parent` to which the `query` is isomorphic---is the length of `search.isomorphisms`

```jldoctest find
nb_locations(search) # = length(search.isomorphisms)
# output
24
```

Element `i_loc` of `search.isomorphisms`, `search.isomorphisms[i_loc]`, is a vector of isomorphisms that share the same subset of atoms in the `parent`, each of which correspond to a different orientation of the `query` overlaying that `parent` substructure. The function `nb_ori_at_loc` outputs a vector whose element `i_loc` is the number of overlay orientations at that location.

```jldoctest find; output=false
nb_ori_at_loc(search)  # 24-element Vector{Int64}: [4, 4, 4, ..., 4]
# output
24-element Vector{Int64}:
 4
 4
 4
 4
 4
 4
 4
 4
 4
 4
 ⋮
 4
 4
 4
 4
 4
 4
 4
 4
 4
```

The individual isomorphisms `isom = search.isomorphisms[i_loc][i_ori]` for a specific location `i_loc` and orientation `i_ori` indicate the correspondence from the `query` to the `parent` struture.
If atom `q` of the `query` maps to atom `p` of the parent, then `isom[q] == p`.

The total number of isomorphisms is given by `nb_isomorphisms(search)`.

```jldoctest find
nb_isomorphisms(search) # = sum(nb_ori_at_loc(search))
# output
96
```

N.b. to generate a `Crystal` containing only the substructures of the `parent` which are isomorphic to the `query`, use:

```jldoctest find
isomorphic_substructures(search)
# output
Name: IRMOF-1.cif
Bravais unit cell of a crystal.
	Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
	Unit cell dimensions a = 25.832000 Å. b = 25.832000 Å, c = 25.832000 Å
	Volume of unit cell: 17237.492730 Å³

	# atoms = 240
	# charges = 0
	chemical formula: Dict(:H => 2, :C => 3)
	space Group: P1
	symmetry Operations:
		'x, y, z'
```

## Molecular and Graph Symmetry

Due to the representation of molecules as graphs, `PoreMatMod.jl` may find more subgraph matches than you may at first expect. For example, searching for CH$_3$ in C$_2$H$_6$ yields not two matches, but six matches because there are three ways we can overlay CH$_3$ on each CH$_4$ group of the ethane.
Sometimes, these subgraphs correspond to distinct spatial isomers or orientations.
We advise to define the `query` using the most minimal structure that matches the targeted `parent` substructure.

An example is searching for [BDC.xyz](../../../assets/find/BDC.xyz) `query` in IRMOF-1 `parent` instead of the more minimal *p*-phenylene fragment.
Thanks to the two carboxyl groups, the total number of isomorphisms is multiplied by a factor of 4, due to the graph-equivalence of the oxygen atoms in each group.  
The number of _locations_ at which the isomorphisms are found, however, is unchanged.

```jldoctest find
query = moiety("BDC.xyz")
search = query ∈ parent
nb_isomorphisms(search) 
# output
384
```
```jldoctest find
nb_locations(search) 
# output
24
```

## Documentation for functions

```@docs
Search
substructure_search
nb_ori_at_loc
nb_isomorphisms
nb_locations
isomorphic_substructures
```
