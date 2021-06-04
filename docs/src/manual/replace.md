```@meta
DocTestSetup = quote
    using MOFun
    xtal = Crystal("IRMOF-1.cif")
    infer_bonds!(xtal, true)
end
```

# Find/Replace Operations

In the previous sections, we saw how we can represent structures using the
[`Xtals`](https://github.com/SimonEnsemble/Xtals.jl) framework, and how we can
identify substructure matches via Ullmann's algorithm for subgraph isomorphism.
Next, we will see how to use this to effect a substructure replacement.

## Orthogonal Procrustes

The question of how best to align crystal fragments in space is addressed by the
["Orthogonal Procrustes"](https://simonensemble.github.io/2018-10/orthogonal-procrustes.html)
problem. If two point sets $A$ and $B$ have a known 1-to-1 correspondence mapping
each point in $A$ to a unique point in $B$, and $A$ is a "noisy" rotated copy of
$B$, we wish to find $R$, the rotation matrix that transforms $A$ to reconstruct
$B$ with the minimum error. This is accomplished via singular value decomposition
of $ABᵀ$. In the linked example $A$ and $B$ are two-dimensional for the purpose
of visualization, but the solution to the problem is the same for point sets in
any $n$-dimensional space. In `MOFun` the solution is applied three-dimensionally,
using the subgraph isomorphism identified by Ullmann's algorithm as the 1-to-1 correspondence.

## The Find-Replace Algorithm

With Ullmann's Algorithm (UA) and Orthogonal Procrustes (OP) in mind, we can tackle the
problem of replacing the BDC linkers of IRMOF-1 with the 2-acetylamido derivative.
Having performed a search for 2-!-*p*-phenylene and chosen a single isomorphism from
the results, we have the 1-to-1 correspondence we need to perform OP to align our
search moiety's notion of Cartesian space with the unit cell of our crystal.  This
correspondence will be termed $s→p$. We also need OP to align the replacement moiety,
2-acetylamido-*p*-phenylene, with the search moiety, so we perform a second UA.

The search moiety in this case is the original search moiety minus the subset of
atoms tagged with `!`. For 2-!-*p*-phenylene, this new search term is
1,2,4-dehydro-benzene, by deletion of the `:H!` at the 2 position of *p*-phenylene.
Performing the second UA (seeking 1,2,4-dehydro-benzene in 2-acetylamido-*p*-phenylene)
affords a second 1-to-1 correspondence, $r→s′$.

OP with $s→p$ and $r→s′$ give two rotation matrices, which transform the replacement
moiety by sequential multiplication to align it with the coordinates of the parent
crystal such that it can replace the linker at the chosen location. The transformed
replacement moiety is added to the crystal, overwriting the original location's atoms.

`MOFun` does all of this in one line of user code. Bonds, including any across periodic
boundaries, are preserved, the unit cell's dimensions are maintained, and the new
structure can be saved to disk for use in simulations.

## Inputs

To guide replacement, the `.xyz` file input for the search moiety must be
altered.  Simply adding `!` after the atomic symbol tags that atom for replacement.

In the example of finding and replacing the 2 position H of the *p*-phenylene moieties
of IRMOF-1, the input should have one H atom tagged, like so:

```
H!         1.06706        0.70670        1.48683
```

Load the new file in like before.

```jldoctest replace_md; output=false
s_moty = moiety("2-!-p-phenylene")
search = s_moty in xtal
# output
2-!-p-phenylene ∈ IRMOF-1.cif
96 hits in 24 locations.
```

The `!` tag does not affect the outcome of [`substructure_search`](@ref).

```jldoctest
s1 = substructure_search(moiety("p-phenylene"), xtal) 
s2 = substructure_search(moiety("2-!-p-phenylene"), xtal)

nb_isomorphisms(s1) == nb_isomorphisms(s2) &&
nb_locations(s1) == nb_locations(s2) &&
isequal(nb_configs_at_loc(s1), nb_configs_at_loc(s2))
# output
true
```

A final additional file input is required: a `.xyz` file containing the replacement
moiety.  This structure must contain a subgraph isomorphic to the non-tagged atoms
of the search moiety.

To replace *p*-phenylene moieties with 2-acetylamido-*p*-phenylene moieties, provide
the [appropriate data](assets/2-acetylamido-p-phenylene) in [`rc[:paths][:moieties]`](manual/inputs)
and load the replacement moiety:

```jldoctest replace_md
r_moty = moiety("2-acetylamido-p-phenylene")
# output
Name: 2-acetylamido-p-phenylene
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

With all three file inputs loaded (parent crystal IRMOF-1 as `xtal`, search moiety
2-!-*p*-phenylene as`s_moty`, and replacement moiety 2-acetylamido-*p*-phenylene as
`r_moty`) and a substructure search `search` performed, replacements can be made.

`MOFun.jl` has several replacement modes, one of which must be specified.

### Random orientation at each location

Random configurations will be chosen for each location in `search.results`, so
that each occurrence of the search moiety in the parent crystal is replaced.

```jldoctest replace_md
substructure_replace(search, r_moty, rand_all=true)
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

### Random orientation at n random locations

The parent crystal's search moieties will be replaced using random configurations
at each of $n$ locations.

```jldoctest replace_md
substructure_replace(search, r_moty, nb_loc=4)
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

### Random orientation at specific locations

The parent crystal is modified using random configurations at a list of specified
locations.

```jldoctest replace_md
substructure_replace(search, r_moty, loc=[13, 20])
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
# output

```

### Specific orientations at specific locations

Specific replacements are made.

```jldoctest replace_md
substructure_replace(search, r_moty, loc=[13, 20], ori=[1, 1])
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

## Documentation

```@docs
replace
```
