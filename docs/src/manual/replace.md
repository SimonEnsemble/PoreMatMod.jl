```@meta
DocTestSetup = quote
    using PoreMatMod
    parent = Crystal("IRMOF-1.cif")
    infer_bonds!(parent, true)
end
```

# Find/Replace Operations

Next, we learn to effect a substructure replacement. To learn by example, we tackle the problem of replacing the BDC linkers of IRMOF-1 with the 2-acetylamido-BDC linker.

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
isequal(nb_ori_at_loc(s1), nb_ori_at_loc(s2))
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

## Orthogonal Procrustes

The question of how best to align crystal fragments in space is addressed by the ["Orthogonal Procrustes"](https://simonensemble.github.io/2018-10/orthogonal-procrustes.html) problem. 
If two point sets $A$ and $B$ have a known 1-to-1 correspondence mapping each point in $A$ to a unique point in $B$, and $A$ is a "noisy" rotated copy of $B$, we wish to find $R$, the rotation matrix that transforms $A$ to reconstruct $B$ with the minimum error. 
This is accomplished via singular value decomposition of $ABᵀ$. 
In the linked example $A$ and $B$ are two-dimensional for the purpose of visualization, but the solution to the problem is the same for point sets in any $n$-dimensional space. 
In `PoreMatMod.jl` the solution is applied three-dimensionally, using the subgraph isomorphism identified by Ullmann's algorithm as the 1-to-1 correspondence.

## Documentation of functions

```@docs
substructure_replace
```
