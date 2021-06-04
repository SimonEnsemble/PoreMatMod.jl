```@meta
DocTestSetup = quote
    using MOFun
    xtal = Crystal("IRMOF-1.cif")
    infer_bonds!(xtal, true)
    s_moty = moiety("p-phenylene")
end
```

# Substructure Searches

Having seen how to load inputs to `MOFun`, we can take the next step towards
accomplishing our example task of functionalizing IRMOF-1. In order to identify
the *p*-phenylene moiety in IRMOF-1, a variation of Ullmann's algorithm is applied.

## Ullmann's Algorithm

[Ullmann's algorithm for subgraph isomorphism](https://doi.org/10.1145/321921.321925)
is the basis of the substructure search in `MOFun`. The algorithm is a depth-first
search of the permutation tree for all possible one-to-one correspondences between the
nodes of one graph (the search graph) and any subset of nodes of another graph (the
parent graph). The search tree is greatly reduced in size by imposing several constraints
on possible node-to-node correspondences. At each branch of the search tree, additional 
pruning further reduces the search space by comparing the immediate neighborhoods of 
potentially-corresponding nodes.

`MOFun` augments Ullmann's algorithm to include the requirement that potentially-
corresponding nodes be of the same atomic species, as required by the chemical
structure application. Additionally, the number of initial potential matches is
reduced by further examination of each atom's local bonding neighborhood.  Bonds
across the unit cell boundaries of periodic materials are handled innately by
the graph representation of structures in `Xtals`.

## Searching with `MOFun`

With a parent crystal and search moiety loaded, execute a search:

```jldoctest
search = substructure_search(s_moty, xtal)
# output
p-phenylene ∈ IRMOF-1.cif
96 hits in 24 locations.
```

The `∈` (`in`) infix operator will also perform the search:

```jldoctest find
search = s_moty ∈ xtal
# output
p-phenylene ∈ IRMOF-1.cif
96 hits in 24 locations.
```

This returns a `Search` object.  Its `query` attribute stores the `Crystal` data
for the search and parent structures, and its `results` attribute is a
`GroupedDataFrame` listing the subgraph isomorphisms grouped by location.

```jldoctest find
search.query
# output
p-phenylene ∈ IRMOF-1.cif
```
```jldoctest find
typeof(search.query.parent)
# output
Crystal
```
```jldoctest find
search.query.s_moty.name
# output
"p-phenylene"
```
```jldoctest find
search.results
# output
GroupedDataFrame with 24 groups based on key: location
First Group (4 rows): location = 1
│ Row │ location │ isomorphism                                        │
│     │ Int64    │ Array{Int64,N} where N                             │
├─────┼──────────┼────────────────────────────────────────────────────┤
│ 1   │ 1        │ [233, 306, 318, 245, 185, 197, 414, 329, 402, 341] │
│ 2   │ 1        │ [245, 318, 306, 233, 197, 185, 402, 341, 414, 329] │
│ 3   │ 1        │ [306, 233, 245, 318, 185, 197, 341, 402, 329, 414] │
│ 4   │ 1        │ [318, 245, 233, 306, 197, 185, 329, 414, 341, 402] │
⋮
Last Group (4 rows): location = 24
│ Row │ location │ isomorphism                                        │
│     │ Int64    │ Array{Int64,N} where N                             │
├─────┼──────────┼────────────────────────────────────────────────────┤
│ 1   │ 24       │ [288, 311, 326, 301, 228, 229, 422, 384, 407, 397] │
│ 2   │ 24       │ [301, 326, 311, 288, 229, 228, 407, 397, 422, 384] │
│ 3   │ 24       │ [311, 288, 301, 326, 228, 229, 397, 407, 384, 422] │
│ 4   │ 24       │ [326, 301, 288, 311, 229, 228, 384, 422, 397, 407] │
```

In the chosen example, the search moiety (*p*-phenylene) occurs 24 times in the
provided structure of the parent crystal (IRMOF-1), with 4 symmetry-equivalent search
hits at each site, for a total of 96 subgraph isomorphisms.

```jldoctest find
nb_locations(search) 
# output
24
```jldoctest find
nb_configs_at_loc(search) 
# output
[4, 4, 4, ..., 4]
```
```jldoctest find
nb_isomorphisms(search) 
# output
96
```

## Symmetry and Rotational Freedom

Due to the representation of molecules as graphs, independently rotatable groups
may yield "extra" search results corresponding to different rotamers.  In some
applications this may be advantageous, but in most cases it is advisable to search
using the most minimal structure which uniquely matches the targeted parent moiety.

An example is searching for [*p*-terephthalate](../../../assets/p-terephthalate.xyz)
in IRMOF-1 instead of the more minimal *p*-phenylene.  Thanks to the two
independently rotatable carboxyl groups, the total number of isomorphisms is
multiplied by a factor of 4.  The number of locations at which the isomorphisms
are found is unchanged.

```jldoctest find
s_moty = moiety("p-terephthalate")
search = s_moty ∈ xtal
nb_isomorphisms(search) 
# output
384
```
```jldoctest find
nb_locations(search) 
# output
24
```

## Documentation

```@docs
Search
Query
substructure_search
nb_configs_at_loc
nb_isomorphisms
nb_locations
```
