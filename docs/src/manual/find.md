```@meta
DocTestSetup = quote
    using PoreMatMod
    parent = Crystal("IRMOF-1.cif")
    infer_bonds!(parent, true)
    query = moiety("p-phenylene.xyz")
end
```

# Subgraph matching (substructure searches)

`PoreMatMod.jl` conducts subgraph matching, i.e. searches for subgraphs of a `parent` graph isomorphic to a `query` graph, using [Ullmann's algorithm for subgraph isomorphisms](https://doi.org/10.1145/321921.321925).

For subgraph matching, both the `parent` crystal structure and `query` fragment are represented by node-labeled (by the chemical species) graphs (nodes = atoms, edges = bonds). For crystals, bonds across the unit cell boundaries of periodic materials are accounted for, allowing us to find subgraph isomorphisms when the fragment is split across a unit cell boundary.

# Substructure Searches: how to

To learn by example, suppose we wish to search the IRMOF-1 crystal structure for *p*-phenylene fragments.

![find graphic](../../assets/find/s_moty-in-xtal.png)

First, we load the `query` fragment and `parent` crystal structure:

```julia
parent = Crystal("IRMOF-1.cif")
infer_bonds!(parent_xtal, true) # true to infer bonds across the periodic boundary

query = moiety("p-phenylene.xyz")
```

then execute a search for subgraphs of the `parent` that "match" (are isomorphic to) the graph of the `query` fragment:

```jldoctest find
search = substructure_search(query, parent)

# output

p-phenylene.xyz ∈ IRMOF-1.cif
96 hits in 24 locations.
```

!!! note "Syntactic sugar for substructure search"
    
    The `∈` (`\in` then hit `Tab` for this Unicode character) or `in` infix operators will also execute the search:
    
    ```jldoctest find; output=false
    search = query ∈ parent
    # or
    search = query in parent
    # or
    search = substructure_search(query, parent)
    
    # output
    
    p-phenylene.xyz ∈ IRMOF-1.cif
    96 hits in 24 locations.
    ```

Both functions `substructure_search` and `∈` return a `Search` object with attributes:

  - `search.query`: the query in the search
  - `search.parent`: the parent in the search
  - `search.isomorphisms`: the result of a search---a nested vector giving the query-to-parent correpondence dictionaries.

```jldoctest find
search.isomorphisms

# output

24-element Vector{Vector{Dict{Int64, Int64}}}:
 [Dict(5 => 185, 4 => 245, 6 => 197, 7 => 414, 2 => 306, 10 => 341, 9 => 402, 8 => 329, 3 => 318, 1 => 233…), Dict(5 => 197, 4 => 233, 6 => 185, 7 => 402, 2 => 318, 10 => 329, 9 => 414, 8 => 341, 3 => 306, 1 => 245…), Dict(5 => 185, 4 => 318, 6 => 197, 7 => 341, 2 => 233, 10 => 414, 9 => 329, 8 => 402, 3 => 245, 1 => 306…), Dict(5 => 197, 4 => 306, 6 => 185, 7 => 329, 2 => 245, 10 => 402, 9 => 341, 8 => 414, 3 => 233, 1 => 318…)]
 [Dict(5 => 186, 4 => 246, 6 => 198, 7 => 413, 2 => 305, 10 => 342, 9 => 401, 8 => 330, 3 => 317, 1 => 234…), Dict(5 => 198, 4 => 234, 6 => 186, 7 => 401, 2 => 317, 10 => 330, 9 => 413, 8 => 342, 3 => 305, 1 => 246…), Dict(5 => 186, 4 => 317, 6 => 198, 7 => 342, 2 => 234, 10 => 413, 9 => 330, 8 => 401, 3 => 246, 1 => 305…), Dict(5 => 198, 4 => 305, 6 => 186, 7 => 330, 2 => 246, 10 => 401, 9 => 342, 8 => 413, 3 => 234, 1 => 317…)]
 [Dict(5 => 187, 4 => 248, 6 => 200, 7 => 415, 2 => 308, 10 => 344, 9 => 404, 8 => 331, 3 => 319, 1 => 235…), Dict(5 => 200, 4 => 235, 6 => 187, 7 => 404, 2 => 319, 10 => 331, 9 => 415, 8 => 344, 3 => 308, 1 => 248…), Dict(5 => 187, 4 => 319, 6 => 200, 7 => 344, 2 => 235, 10 => 415, 9 => 331, 8 => 404, 3 => 248, 1 => 308…), Dict(5 => 200, 4 => 308, 6 => 187, 7 => 331, 2 => 248, 10 => 404, 9 => 344, 8 => 415, 3 => 235, 1 => 319…)]
 [Dict(5 => 188, 4 => 247, 6 => 199, 7 => 416, 2 => 307, 10 => 343, 9 => 403, 8 => 332, 3 => 320, 1 => 236…), Dict(5 => 199, 4 => 236, 6 => 188, 7 => 403, 2 => 320, 10 => 332, 9 => 416, 8 => 343, 3 => 307, 1 => 247…), Dict(5 => 188, 4 => 320, 6 => 199, 7 => 343, 2 => 236, 10 => 416, 9 => 332, 8 => 403, 3 => 247, 1 => 307…), Dict(5 => 199, 4 => 307, 6 => 188, 7 => 332, 2 => 247, 10 => 403, 9 => 343, 8 => 416, 3 => 236, 1 => 320…)]
 [Dict(5 => 189, 4 => 255, 6 => 207, 7 => 376, 2 => 262, 10 => 351, 9 => 358, 8 => 333, 3 => 280, 1 => 237…), Dict(5 => 207, 4 => 237, 6 => 189, 7 => 358, 2 => 280, 10 => 333, 9 => 376, 8 => 351, 3 => 262, 1 => 255…), Dict(5 => 189, 4 => 280, 6 => 207, 7 => 351, 2 => 237, 10 => 376, 9 => 333, 8 => 358, 3 => 255, 1 => 262…), Dict(5 => 207, 4 => 262, 6 => 189, 7 => 333, 2 => 255, 10 => 358, 9 => 351, 8 => 376, 3 => 237, 1 => 280…)]
 [Dict(5 => 190, 4 => 256, 6 => 208, 7 => 375, 2 => 261, 10 => 352, 9 => 357, 8 => 334, 3 => 279, 1 => 238…), Dict(5 => 208, 4 => 238, 6 => 190, 7 => 357, 2 => 279, 10 => 334, 9 => 375, 8 => 352, 3 => 261, 1 => 256…), Dict(5 => 190, 4 => 279, 6 => 208, 7 => 352, 2 => 238, 10 => 375, 9 => 334, 8 => 357, 3 => 256, 1 => 261…), Dict(5 => 208, 4 => 261, 6 => 190, 7 => 334, 2 => 256, 10 => 357, 9 => 352, 8 => 375, 3 => 238, 1 => 279…)]
 [Dict(5 => 191, 4 => 254, 6 => 206, 7 => 373, 2 => 264, 10 => 350, 9 => 360, 8 => 335, 3 => 277, 1 => 239…), Dict(5 => 206, 4 => 239, 6 => 191, 7 => 360, 2 => 277, 10 => 335, 9 => 373, 8 => 350, 3 => 264, 1 => 254…), Dict(5 => 191, 4 => 277, 6 => 206, 7 => 350, 2 => 239, 10 => 373, 9 => 335, 8 => 360, 3 => 254, 1 => 264…), Dict(5 => 206, 4 => 264, 6 => 191, 7 => 335, 2 => 254, 10 => 360, 9 => 350, 8 => 373, 3 => 239, 1 => 277…)]
 [Dict(5 => 192, 4 => 253, 6 => 205, 7 => 374, 2 => 263, 10 => 349, 9 => 359, 8 => 336, 3 => 278, 1 => 240…), Dict(5 => 205, 4 => 240, 6 => 192, 7 => 359, 2 => 278, 10 => 336, 9 => 374, 8 => 349, 3 => 263, 1 => 253…), Dict(5 => 192, 4 => 278, 6 => 205, 7 => 349, 2 => 240, 10 => 374, 9 => 336, 8 => 359, 3 => 253, 1 => 263…), Dict(5 => 205, 4 => 263, 6 => 192, 7 => 336, 2 => 253, 10 => 359, 9 => 349, 8 => 374, 3 => 240, 1 => 278…)]
 [Dict(5 => 193, 4 => 252, 6 => 204, 7 => 395, 2 => 290, 10 => 348, 9 => 386, 8 => 337, 3 => 299, 1 => 241…), Dict(5 => 204, 4 => 241, 6 => 193, 7 => 386, 2 => 299, 10 => 337, 9 => 395, 8 => 348, 3 => 290, 1 => 252…), Dict(5 => 193, 4 => 299, 6 => 204, 7 => 348, 2 => 241, 10 => 395, 9 => 337, 8 => 386, 3 => 252, 1 => 290…), Dict(5 => 204, 4 => 290, 6 => 193, 7 => 337, 2 => 252, 10 => 386, 9 => 348, 8 => 395, 3 => 241, 1 => 299…)]
 [Dict(5 => 194, 4 => 251, 6 => 203, 7 => 396, 2 => 289, 10 => 347, 9 => 385, 8 => 338, 3 => 300, 1 => 242…), Dict(5 => 203, 4 => 242, 6 => 194, 7 => 385, 2 => 300, 10 => 338, 9 => 396, 8 => 347, 3 => 289, 1 => 251…), Dict(5 => 194, 4 => 300, 6 => 203, 7 => 347, 2 => 242, 10 => 396, 9 => 338, 8 => 385, 3 => 251, 1 => 289…), Dict(5 => 203, 4 => 289, 6 => 194, 7 => 338, 2 => 251, 10 => 385, 9 => 347, 8 => 396, 3 => 242, 1 => 300…)]
 ⋮
 [Dict(5 => 212, 4 => 271, 6 => 219, 7 => 392, 2 => 283, 10 => 367, 9 => 379, 8 => 356, 3 => 296, 1 => 260…), Dict(5 => 219, 4 => 260, 6 => 212, 7 => 379, 2 => 296, 10 => 356, 9 => 392, 8 => 367, 3 => 283, 1 => 271…), Dict(5 => 212, 4 => 296, 6 => 219, 7 => 367, 2 => 260, 10 => 392, 9 => 356, 8 => 379, 3 => 271, 1 => 283…), Dict(5 => 219, 4 => 283, 6 => 212, 7 => 356, 2 => 271, 10 => 379, 9 => 367, 8 => 392, 3 => 260, 1 => 296…)]
 [Dict(5 => 213, 4 => 276, 6 => 224, 7 => 419, 2 => 314, 10 => 372, 9 => 410, 8 => 361, 3 => 323, 1 => 265…), Dict(5 => 224, 4 => 265, 6 => 213, 7 => 410, 2 => 323, 10 => 361, 9 => 419, 8 => 372, 3 => 314, 1 => 276…), Dict(5 => 213, 4 => 323, 6 => 224, 7 => 372, 2 => 265, 10 => 419, 9 => 361, 8 => 410, 3 => 276, 1 => 314…), Dict(5 => 224, 4 => 314, 6 => 213, 7 => 361, 2 => 276, 10 => 410, 9 => 372, 8 => 419, 3 => 265, 1 => 323…)]
 [Dict(5 => 214, 4 => 275, 6 => 223, 7 => 420, 2 => 313, 10 => 371, 9 => 409, 8 => 362, 3 => 324, 1 => 266…), Dict(5 => 223, 4 => 266, 6 => 214, 7 => 409, 2 => 324, 10 => 362, 9 => 420, 8 => 371, 3 => 313, 1 => 275…), Dict(5 => 214, 4 => 324, 6 => 223, 7 => 371, 2 => 266, 10 => 420, 9 => 362, 8 => 409, 3 => 275, 1 => 313…), Dict(5 => 223, 4 => 313, 6 => 214, 7 => 362, 2 => 275, 10 => 409, 9 => 371, 8 => 420, 3 => 266, 1 => 324…)]
 [Dict(5 => 215, 4 => 273, 6 => 221, 7 => 418, 2 => 316, 10 => 369, 9 => 412, 8 => 363, 3 => 322, 1 => 267…), Dict(5 => 221, 4 => 267, 6 => 215, 7 => 412, 2 => 322, 10 => 363, 9 => 418, 8 => 369, 3 => 316, 1 => 273…), Dict(5 => 215, 4 => 322, 6 => 221, 7 => 369, 2 => 267, 10 => 418, 9 => 363, 8 => 412, 3 => 273, 1 => 316…), Dict(5 => 221, 4 => 316, 6 => 215, 7 => 363, 2 => 273, 10 => 412, 9 => 369, 8 => 418, 3 => 267, 1 => 322…)]
 [Dict(5 => 216, 4 => 274, 6 => 222, 7 => 417, 2 => 315, 10 => 370, 9 => 411, 8 => 364, 3 => 321, 1 => 268…), Dict(5 => 222, 4 => 268, 6 => 216, 7 => 411, 2 => 321, 10 => 364, 9 => 417, 8 => 370, 3 => 315, 1 => 274…), Dict(5 => 216, 4 => 321, 6 => 222, 7 => 370, 2 => 268, 10 => 417, 9 => 364, 8 => 411, 3 => 274, 1 => 315…), Dict(5 => 222, 4 => 315, 6 => 216, 7 => 364, 2 => 274, 10 => 411, 9 => 370, 8 => 417, 3 => 268, 1 => 321…)]
 [Dict(5 => 225, 4 => 303, 6 => 231, 7 => 424, 2 => 310, 10 => 399, 9 => 406, 8 => 381, 3 => 328, 1 => 285…), Dict(5 => 231, 4 => 285, 6 => 225, 7 => 406, 2 => 328, 10 => 381, 9 => 424, 8 => 399, 3 => 310, 1 => 303…), Dict(5 => 225, 4 => 328, 6 => 231, 7 => 399, 2 => 285, 10 => 424, 9 => 381, 8 => 406, 3 => 303, 1 => 310…), Dict(5 => 231, 4 => 310, 6 => 225, 7 => 381, 2 => 303, 10 => 406, 9 => 399, 8 => 424, 3 => 285, 1 => 328…)]
 [Dict(5 => 226, 4 => 304, 6 => 232, 7 => 423, 2 => 309, 10 => 400, 9 => 405, 8 => 382, 3 => 327, 1 => 286…), Dict(5 => 232, 4 => 286, 6 => 226, 7 => 405, 2 => 327, 10 => 382, 9 => 423, 8 => 400, 3 => 309, 1 => 304…), Dict(5 => 226, 4 => 327, 6 => 232, 7 => 400, 2 => 286, 10 => 423, 9 => 382, 8 => 405, 3 => 304, 1 => 309…), Dict(5 => 232, 4 => 309, 6 => 226, 7 => 382, 2 => 304, 10 => 405, 9 => 400, 8 => 423, 3 => 286, 1 => 327…)]
 [Dict(5 => 227, 4 => 302, 6 => 230, 7 => 421, 2 => 312, 10 => 398, 9 => 408, 8 => 383, 3 => 325, 1 => 287…), Dict(5 => 230, 4 => 287, 6 => 227, 7 => 408, 2 => 325, 10 => 383, 9 => 421, 8 => 398, 3 => 312, 1 => 302…), Dict(5 => 227, 4 => 325, 6 => 230, 7 => 398, 2 => 287, 10 => 421, 9 => 383, 8 => 408, 3 => 302, 1 => 312…), Dict(5 => 230, 4 => 312, 6 => 227, 7 => 383, 2 => 302, 10 => 408, 9 => 398, 8 => 421, 3 => 287, 1 => 325…)]
 [Dict(5 => 228, 4 => 301, 6 => 229, 7 => 422, 2 => 311, 10 => 397, 9 => 407, 8 => 384, 3 => 326, 1 => 288…), Dict(5 => 229, 4 => 288, 6 => 228, 7 => 407, 2 => 326, 10 => 384, 9 => 422, 8 => 397, 3 => 311, 1 => 301…), Dict(5 => 228, 4 => 326, 6 => 229, 7 => 397, 2 => 288, 10 => 422, 9 => 384, 8 => 407, 3 => 301, 1 => 311…), Dict(5 => 229, 4 => 311, 6 => 228, 7 => 384, 2 => 301, 10 => 407, 9 => 397, 8 => 422, 3 => 288, 1 => 326…)]
```

In this example, the `query` fragment (*p*-phenylene) occurs in 24 different location in the `parent` crystal structure, with 4 symmetry-equivalent isomorphisms at each location, for a total of 96 subgraph isomorphisms.

The number of locations---the number of unique substructures of the `parent` to which the `query` is isomorphic---is the length of `search.isomorphisms`.

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

Each individual isomorphism `isom = search.isomorphisms[i_loc][i_ori]` for a specific location `i_loc` and orientation `i_ori` indicates the correspondence from the `query` to the `parent` struture: if atom `q` of the `query` maps to atom `p` of the `parent`, then `isom[q] == p`.

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

Crystal(C₁₄₄H₉₆, periodic = TTT):
    bounding_box      : [  25.832        0        0;
                         1.58175e-15   25.832        0;
                         1.58175e-15 1.58175e-15   25.832]u"Å"
```

## Stereochemistry and Isomorphism

The node-labeled graph representation of a molecule/crystal structure is invariant with respect to stereochemistry.
In other words, every rotational/conformational state and stereoisomer of a structure share the same graph representation.
What this means is that `PoreMatMod.jl` may find more subgraph matches than you may first expect.

*Example 1*: Suppose we search for a carboxylate with beta hydrogen in acrylate.

![symmetry viz](../../assets/find/symmetry.png)

There is clearly only one substructure of acrylate that matches the query.
However, there are two subgraph isomorphisms, because swapping the oxygen atoms in the point cloud representation results in the same graph representation.
The above image gives a closer look at how these degenerate representations translate to multiple isomorphisms for a single occurence of a fragment in a structure.

*Example 2*: Suppose we search the IRMOF-1 `parent` structure for the [BDC.xyz](../../../assets/find/BDC.xyz) linker as the `query` instead of the more minimal *p*-phenylene `query` fragment.
Thanks to the two carboxyl groups, the total number of isomorphisms is multiplied by a factor of 4, due to the 180 degree rotation of these groups having no effect on the graph representation.
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

**Note**: We advise to define the `query` using the most minimal structure that matches the targeted `parent` substructure.

## Documentation for functions

```@docs
Search
substructure_search
nb_ori_at_loc
nb_isomorphisms
nb_locations
isomorphic_substructures
```
