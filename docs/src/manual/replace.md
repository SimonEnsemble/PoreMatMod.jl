```@meta
DocTestSetup = quote
    using PoreMatMod
    parent = Crystal("IRMOF-1.cif")
    infer_bonds!(parent, true)
end
```

# Find/Replace Operations

In the previous sections, we saw how we can represent structures using the [`Xtals`](https://github.com/SimonEnsemble/Xtals.jl) framework, and how we can identify substructure matches via Ullmann's algorithm for subgraph isomorphism.
Next, we will see how to use this to effect a substructure replacement.

## Orthogonal Procrustes

The question of how best to align crystal fragments in space is addressed by the ["Orthogonal Procrustes"](https://simonensemble.github.io/2018-10/orthogonal-procrustes.html) problem. 
If two point sets $A$ and $B$ have a known 1-to-1 correspondence mapping each point in $A$ to a unique point in $B$, and $A$ is a "noisy" rotated copy of $B$, we wish to find $R$, the rotation matrix that transforms $A$ to reconstruct $B$ with the minimum error. 
This is accomplished via singular value decomposition of $ABᵀ$. 
In the linked example $A$ and $B$ are two-dimensional for the purpose of visualization, but the solution to the problem is the same for point sets in any $n$-dimensional space. 
In `PoreMatMod.jl` the solution is applied three-dimensionally, using the subgraph isomorphism identified by Ullmann's algorithm as the 1-to-1 correspondence.

## The Find-Replace Algorithm

With Ullmann's Algorithm (UA) and Orthogonal Procrustes (OP) in mind, we can tackle the problem of replacing the BDC linkers of IRMOF-1 with the 2-acetylamido derivative.

First, we perform a search using a `query` with `!`-tagged atom(s), to direct the replacement stage of the algorithm.
In this case, we should use [2-!-p-phenylene.xyz](../../../assets/replace/2-!-p-phenylene.xyz).
Having performed a search and chosen one or more isomorphisms from the results, we have the 1-to-1 correspondence(s) we need to perform OP to align the coordinate systems of the `parent` and `query`.  
Such a correspondence will be termed $q→p′$. 
We also need to align the `replacement` moiety, 2-acetylamido-*p*-phenylene, with the `query`, so we perform a second UA.

The `query` in this case is the original search moiety minus the subset of atoms tagged with `!`. 
For 2-!-*p*-phenylene, this new search term is 1,2,4-dehydro-benzene, by deletion of the `:H!` at the 2 position of *p*-phenylene.
Performing the second UA (seeking 1,2,4-dehydro-benzene in 2-acetylamido-*p*-phenylene) affords a second 1-to-1 correspondence, $q′→r$, the inverse of which is $r′→q′$.
Chaining these together, we have $r′→p′$, the 1-to-1 correspondence for matching a subset of the `replacement` to a subset of the `parent`.

OP with $r′→q′$ gives a rotation matrix which transforms the `replacement` by multiplication to align it with the coordinates of the parent crystal such that it can replace the linker at the chosen location. 
The transformed replacement moiety is added to the crystal, overwriting the original location's atoms.

`PoreMatMod.jl` does all of this in one line of user code.
Bonds, including any across periodic boundaries, are preserved, the unit cell's dimensions are maintained, and the new structure can be saved to disk for use in simulations.

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
isequal(nb_configs_at_loc(s1), nb_configs_at_loc(s2))
# output
true
```

A final additional file input is required: a `.xyz` file containing the `replacement`.  
This structure must contain a subgraph isomorphic to the non-tagged atoms of the `query`.

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

## Replacement Modes

With all three file inputs loaded (IRMOF-1 as `parent`, 2-!-*p*-phenylene as `query`, and 2-acetylamido-*p*-phenylene as `replacement`) and a `search` performed, replacements can be made.

`PoreMatMod.jl` has several replacement modes, one of which must be specified.

### Optimal orientation at each location

Optimal configurations will be chosen for each location in `search.results`, so that each occurrence of the `query` in the `parent` is replaced with minimal perturbation of the conserved atoms from the parent structure.

```jldoctest replace_md; output=false
substructure_replace(search, replacement)
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

### Optimal orientation at n random locations

The `parent` will be modified using optimal configurations at each of $n$ randomly-selected locations.

```jldoctest replace_md; output=false
substructure_replace(search, replacement, nb_loc=4)
# output
Name: new_xtal
Bravais unit cell of a crystal.
	Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
	Unit cell dimensions a = 25.832000 Å. b = 25.832000 Å, c = 25.832000 Å
	Volume of unit cell: 17237.492730 Å³

	# atoms = 452
	# charges = 0
	chemical formula: Dict(:N => 1, :Zn => 8, :H => 27, :O => 27, :C => 50)
	space Group: P1
	symmetry Operations:
		'x, y, z'
```

### Optimal orientation at specific locations

The `parent` is modified using optimal configurations at a list of specified locations.

```jldoctest replace_md; output=false
substructure_replace(search, replacement, loc=[13, 20])
# output
Name: new_xtal
Bravais unit cell of a crystal.
	Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
	Unit cell dimensions a = 25.832000 Å. b = 25.832000 Å, c = 25.832000 Å
	Volume of unit cell: 17237.492730 Å³

	# atoms = 438
	# charges = 0
	chemical formula: Dict(:N => 1, :Zn => 16, :H => 51, :O => 53, :C => 98)
	space Group: P1
	symmetry Operations:
		'x, y, z'
```

### Specific orientations at specific locations

Specific replacements are made by providing the values of `loc` and `ori`.
If any values in `ori` are zero, the corresponding location will processed with the optimal replacement.

```jldoctest replace_md; output=false
substructure_replace(search, replacement, loc=[13, 20], ori=[1, 0])
# output
Name: new_xtal
Bravais unit cell of a crystal.
	Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
	Unit cell dimensions a = 25.832000 Å. b = 25.832000 Å, c = 25.832000 Å
	Volume of unit cell: 17237.492730 Å³

	# atoms = 438
	# charges = 0
	chemical formula: Dict(:N => 1, :Zn => 16, :H => 51, :O => 53, :C => 98)
	space Group: P1
	symmetry Operations:
		'x, y, z'
```

### Random orientations

By using the `random` keyword argument, the search for optimal alignment can be skipped, with the value of `ori` being selected at random.

```jldoctest replace_md; output=false
substructure_replace(search, replacement, nb_loc=3, random=true)
# output
Name: new_xtal
Bravais unit cell of a crystal.
	Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
	Unit cell dimensions a = 25.832000 Å. b = 25.832000 Å, c = 25.832000 Å
	Volume of unit cell: 17237.492730 Å³

	# atoms = 445
	# charges = 0
	chemical formula: Dict(:N => 3, :Zn => 32, :H => 105, :O => 107, :C => 198)
	space Group: P1
	symmetry Operations:
		'x, y, z'
```

## Simple syntax

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

Note it is still required to specify a replacement style via keyword argument.

## Documentation

```@docs
substructure_replace
```
