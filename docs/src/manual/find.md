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
* `search.query`: the query in the search
* `search.parent`: the parent in the search
* `search.isomorphisms`: the result of a search---a nested vector giving the query-to-parent correpondence dictionaries.

```jldoctest find
search.isomorphisms
# output
asdf
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
