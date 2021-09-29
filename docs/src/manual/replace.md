```@meta
DocTestSetup = quote
    using PoreMatMod
    parent = Crystal("IRMOF-1.cif")
    infer_bonds!(parent, true)
end
```

# Find/Replace Operations

To do a find-and-replace, we perform a search using a `query` with masked atoms.
A masked atom is treated as a normal atom when performing a search between `query` and parent, and is denoted by appending `!` to its chemical species label.
Atoms in the `parent` that are mapped to atoms in the `query` are either deleted or replaced with atoms from the `replacement` fragment (see below).
The masked atoms are needed so that `PoreMatMod` can determine how to align the `replacement` with the targeted subset of the `parent`--masked atoms are ignored when searching for a mapping between the `query` and the `replacement`.

## Inputs

To guide replacement, the `.xyz` file input for the search moiety must be altered.  
Simply adding `!` after the atomic symbol tags that atom for replacement.
The atom property viewer feature in [iRASPA](https://iraspa.org/) is conducive to this task.

In the example of finding and replacing the 2 position H of the *p*-phenylene moieties of IRMOF-1, the input should have one H atom tagged, like so:

```
H!         1.06706        0.70670        1.48683
```

Load the new file in like before.

```jldoctest replace_md; output=false
query = moiety("2-!-p-phenylene.xyz")
search = query in parent
# output
2-!-p-phenylene.xyz ∈ IRMOF-1.cif
96 hits in 24 locations.
```

The `!` tag does not affect the outcome of `substructure_search`.

```jldoctest
s1 = substructure_search(moiety("p-phenylene.xyz"), parent) 
s2 = substructure_search(moiety("2-!-p-phenylene.xyz"), parent)

nb_isomorphisms(s1) == nb_isomorphisms(s2) &&
nb_locations(s1) == nb_locations(s2) &&
isequal(nb_ori_at_loc(s1), nb_ori_at_loc(s2))
# output
true
```

To replace *p*-phenylene moieties with 2-acetylamido-*p*-phenylene moieties, load [2-acetylamido-p-phenylene.xyz](../../../assets/replace/2-acetylamido-p-phenylene.xyz):

```jldoctest replace_md; output=false
replacement = moiety("2-acetylamido-p-phenylene.xyz")
# output
Name: 2-acetylamido-p-phenylene.xyz
Bravais unit cell of a crystal.
	Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
	Unit cell dimensions a = 1.000000 Å. b = 1.000000 Å, c = 1.000000 Å
	Volume of unit cell: 1.000000 Å³

	# atoms = 17
	# charges = 0
	chemical formula: Dict(:N => 1, :H => 7, :O => 1, :C => 8)
	space Group: P1
	symmetry Operations:
		'x, y, z'
```


## Simple Syntax

Generally, it is advisable to perform the `search` and use `substructure_replace`, as multiple replacement tasks can be performed with a single searching step.
The search is usually the slowest step, and it is desirable not to perform it repeatedly.
However, for one-shot find-and-replace operations, the `replace` function syntax from standard `Julia` may be used:

```jldoctest replace_md; output=false
replace(parent, query => replacement)
# output
Name: new_xtal
Bravais unit cell of a crystal.
	Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
	Unit cell dimensions a = 25.832000 Å. b = 25.832000 Å, c = 25.832000 Å
	Volume of unit cell: 17237.492730 Å³

	# atoms = 592
	# charges = 0
	chemical formula: Dict(:N => 3, :Zn => 4, :H => 21, :O => 16, :C => 30)
	space Group: P1
	symmetry Operations:
		'x, y, z'
```

To direct the number, location, and orientation of the replacements made, use the keyword arguments.  These are described in the docstrings and demonstrated in the [examples](../../../examples).

## Documentation of functions

```@docs
substructure_replace
```
