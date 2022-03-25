```@meta
DocTestSetup = quote
    using PoreMatMod
end
```
# Find/Replace Operations

Suppose we wish to conduct the find-and-replace operations illustrated in the figure below, to produce an acetylamido-functionalized IRMOF-1 structure.

![replacement scheme](../../assets/replace/rep_header.png)

#### the `parent` structure
First, we load the `parent` IRMOF-1 structure and infer its bonds.

```jldoctest replace_md; output=false
parent = Crystal("IRMOF-1.cif")
infer_bonds!(parent, true)
# output
true
```

#### the `query` fragment
Next, we define a `query` fragment as a *p*-phenylene moiety.
To guide the replacement, the masked atoms of the `query` fragment must be annotated with `!` in the `.xyz` input file by appending a `!` character at the end of their atomic symbols.
The atom property viewer feature in [iRASPA](https://iraspa.org/) is useful for figuring out which atom(s) to mask.

!!! note
    A masked atom (marked with `!`) in the `query` fragment implies that the corresponding atom of the `parent` crystal structure (i) must be removed [e.g., to make room for replacement with a different functionality] but (ii) does not correspond with an atom on the `replacement` fragment and thus cannot be used in the process of aligning the `replacement` fragment onto the `parent` crystal. 

In our example, in `2-!-p-phenylene.xyz` input file describing our *p*-phenylene `query` fragment, one H atom is masked (see figure above):

```
10

C         -1.71069        0.96969       -0.46280
C         -0.48337        1.30874        0.11690
C         -2.33707       -0.23371       -0.12103
C          0.11757        0.44439        1.03836
C         -0.50881       -0.75900        1.38013
C         -1.73613       -1.09805        0.80043
H!         1.06706        0.70670        1.48683
H          0.00122        2.23972       -0.14750
H         -3.28655       -0.49601       -0.56950
H         -2.22071       -2.02904        1.06484
```

We then read the input file for the `query` fragment.

```jldoctest replace_md; output=false
query = moiety("2-!-p-phenylene.xyz")
# output
Name: 2-!-p-phenylene.xyz
Bravais unit cell of a crystal.
	Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
	Unit cell dimensions a = 1.000000 Å. b = 1.000000 Å, c = 1.000000 Å
	Volume of unit cell: 1.000000 Å³

	# atoms = 10
	# charges = 0
	chemical formula: C₆H!H₃
	space Group: P1
	symmetry Operations:
		'x, y, z'
	bonding graph:
		# vertices = 10
		# edges = 10
```

#### the `replacement` fragment

Next, we read in the acetylamido-functionalized version of the `query` fragment, [2-acetylamido-p-phenylene.xyz](../../../assets/replace/2-acetylamido-p-phenylene.xyz), as the `replacement` fragment:

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
	chemical formula: C₈H₇NO
	space Group: P1
	symmetry Operations:
		'x, y, z'
	bonding graph:
		# vertices = 17
		# edges = 17
```

#### the find step

We search for subgraphs of the `parent` structure that match the `query` fragment.
Note the `!` tags are ignored during the `substructure_search`.

```jldoctest replace_md; output=false
search = query in parent # equivalent to substructure_search(query, parent)
# output
2-!-p-phenylene.xyz ∈ IRMOF-1.cif
96 hits in 24 locations.
```

#### the replace step
The code below will, at each location in the `parent` where a substructure matched the `query` fragment, choose a random orientation (corresponding to an overlay of the `query` with the substructure), align and install the replacement fragment, then remove the original substructure, giving the `child` structure shown in the figure above.

```jldoctest replace_md; output=false
child = substructure_replace(search, replacement)
# output
Name: new_xtal
Bravais unit cell of a crystal.
	Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
	Unit cell dimensions a = 25.832000 Å. b = 25.832000 Å, c = 25.832000 Å
	Volume of unit cell: 17237.492730 Å³

	# atoms = 592
	# charges = 0
	chemical formula: C₂₄₀H₁₆₈N₂₄O₁₂₈Zn₃₂
	space Group: P1
	symmetry Operations:
		'x, y, z'
	bonding graph:
		# vertices = 592
		# edges = 680
```

To direct the number, location, and orientation of the replacements, use the keyword arguments for [`substructure_replace`](@ref). Particularly, the location `loc` and orientation `ori` keyword arguments specify a particular isomorphism to use (in reference to `search.isomorphisms`) when conducting a replacement operation. The figure below illustrates.

![loc/ori example](../../assets/replace/loc_ori_example.png)

For more details, see the [search docs](../../find) and the [replacement modes example](../../../examples/replacement_modes.html).

### quick find-and-replace syntax

For one-shot find-and-replace operations, the `replace` function may be used:

```jldoctest replace_md; output=false
child = replace(parent, query => replacement)
# output
Name: new_xtal
Bravais unit cell of a crystal.
	Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
	Unit cell dimensions a = 25.832000 Å. b = 25.832000 Å, c = 25.832000 Å
	Volume of unit cell: 17237.492730 Å³

	# atoms = 592
	# charges = 0
	chemical formula: C₂₄₀H₁₆₈N₂₄O₁₂₈Zn₃₂
	space Group: P1
	symmetry Operations:
		'x, y, z'
	bonding graph:
		# vertices = 592
		# edges = 680
```

!!! note
    Generally, it is advisable to perform the search using `substructure_replace` then pass it to `replace`, as multiple replacement tasks can then be performed on the basis of the search step as opposed to repeating it for each replacement. The search is usually the slowest step, and it is desirable not to perform it repeatedly.


## Documentation of functions

```@docs
substructure_replace
replace
```
