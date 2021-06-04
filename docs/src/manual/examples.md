```@meta
DocTestSetup = quote
    using MOFun
end
```

## Examples

Each example case below gives a description, illustration, and code snippet,
and is linked to a more detailed Pluto notebook tutorial.

### Generate hypothetical structures

[![example 1](../../../assets/IRMOF1example.png)
Pluto Notebook](../../../examples/IRMOF1.jl)

Create novel derivatives of a `Crystal`, search for one of its substructures,
and replace with a derivatized moiety.

Example: *ortho* substitution with an acetylamido group at one quarter of the
*p*-phenylene moieties in IRMOF-1.

```jldoctest
xtal = Crystal("IRMOF-1.cif")
infer_bonds!(xtal, true)
s_moty = moiety("2-!-p-phenylene")
r_moty = moiety("2-acetylamido-p-phenylene")
search = s_moty ∈ xtal
new_xtal = substructure_replace(search, r_moty, nb_loc=Int(nb_locations(search)/4))
# output
Name: new_xtal
Bravais unit cell of a crystal.
	Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
	Unit cell dimensions a = 25.832000 Å. b = 25.832000 Å, c = 25.832000 Å
	Volume of unit cell: 17237.492730 Å³

	# atoms = 466
	# charges = 0
	chemical formula: Dict(:N => 3, :Zn => 16, :H => 57, :O => 55, :C => 102)
	space Group: P1
	symmetry Operations:
		'x, y, z'
```

### Insert missing hydrogens

[![example 4](../../../assets/missingHexample.png)
Pluto Notebook](../../../examples/missingH.jl)

To correct defects like missing atoms, use the affected substructure as the search
moiety and a manually corrected copy as the replacement moiety.

Example: Insert missing H atoms in IRMOF-1


```jldoctest
xtal = Crystal("IRMOF-1_noH.cif")
infer_bonds!(xtal, true)
s_moty = moiety("p-phenylene_noH")
r_moty = moiety("p-phenylene")
repaired_xtal = substructure_replace(s_moty ∈ xtal, r_moty, rand_all=true)
# output
Name: new_xtal
Bravais unit cell of a crystal.
	Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
	Unit cell dimensions a = 25.832000 Å. b = 25.832000 Å, c = 25.832000 Å
	Volume of unit cell: 17237.492730 Å³

	# atoms = 424
	# charges = 0
	chemical formula: Dict(:Zn => 4, :H => 12, :O => 13, :C => 24)
	space Group: P1
	symmetry Operations:
		'x, y, z'
```

### Repair Disorder and Remove Adsorbates

[![example 5](../../../assets/landingpageexample.png)
Pluto Notebook](../../../examples/landingpage.jl)

The example on the [landing page](../../../index.md): repair, activate, and functionalize.

Note the use of the `(s_moty => r_moty) in xtal` syntactic sugar.

```jldoctest
using MOFun
xtal = Crystal("EMEHUB_C2H2.cif", remove_duplicates=true, check_overlap=false)
infer_bonds!(xtal, true)
repaired = (moiety("disordered_ligand!") => moiety("4-pyridyl")) ∈ xtal
active = substructure_replace(
    substructure_search(moiety("acetylene"), repaired, exact=true), 
    moiety(nothing), rand_all=true)
# output
┌ Info: Crystal EMEHUB_C2H2.cif has I 4/m m m space group. I am converting it to P1 symmetry.
└         To afrain from this, pass `convert_to_p1=false` to the `Crystal` constructor.
┌ Warning: carbon atom 1 in EMEHUB_C2H2.cif is bonded to more than four atoms!
└ @ Xtals ~/.julia/dev/Xtals/src/bonds.jl:403
┌ Warning: carbon atom 6 in disordered_ligand! is bonded to more than four atoms!
└ @ Xtals ~/.julia/dev/Xtals/src/bonds.jl:403
┌ Warning: carbon atom 1 in disordered_ligand! is bonded to more than four atoms!
└ @ Xtals ~/.julia/dev/Xtals/src/bonds.jl:403
Name: new_xtal
Bravais unit cell of a crystal.
	Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
	Unit cell dimensions a = 13.715000 Å. b = 13.715000 Å, c = 7.952260 Å
	Volume of unit cell: 1495.829848 Å³

	# atoms = 104
	# charges = 0
	chemical formula: Dict(:N => 4, :F => 6, :H => 16, :Cu => 1, :Si => 1, :C => 24)
	space Group: P1
	symmetry Operations:
		'x, y, z'
```
