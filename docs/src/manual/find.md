# Substructure Searches

## Searching

With a parent crystal and search moiety loaded, execute a search:

```
search = substructure_search(s_moty, xtal)
```

The `∈` (`in`) infix operator will also perform the search:

```
search = s_moty ∈ xtal
```

This returns a `Search` object.  Its `query` attribute stores the `Crystal` data
for the search and parent structures, and its `results` attribute is a
`GroupedDataFrame` listing the subgraph isomorphisms grouped by location.

In the chosen example, the search moiety (*p*-phenylene) occurs 24 times in the
parent crystal (IRMOF-1), with 4 symmetry-equivalent search hits at each site,
for a total of 96 subgraph isomorphisms.

```
nb_locations(search) # 24
nb_configs_at_loc(search) # [4, 4, 4, ..., 4]
nb_isomorphisms(search) # 96
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

```
xtal = Crystal("IRMOF-1_clean.cif")
infer_bonds!(xtal, true)
s_moty = moiety("p-terephthalate")
search = s_moty ∈ xtal
nb_isomorphisms(search) # 384
nb_locations(search) # 24
```

## Documentation

```@meta
CurrentModule = MOFun
DocTestSetup = quote
    using MOFun
end
```

```@docs
substructure_search
Search
```
