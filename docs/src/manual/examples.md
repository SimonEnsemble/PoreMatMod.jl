```@meta
DocTestSetup = quote
    using PoreMatMod
end
```

## Examples

### Generate hypothetical structures

Example: *ortho* substitution with an acetylamido group at one quarter of the *p*-phenylene moieties in IRMOF-1.

[Input Files](../../../assets/examples/example1.zip)

```jldoctest; output=false
parent = Crystal("IRMOF-1.cif")
infer_bonds!(parent, true)
query = moiety("2-!-p-phenylene.xyz")
replacement = moiety("2-acetylamido-p-phenylene.xyz")
search = query ∈ parent
new_xtal = substructure_replace(search, replacement, nb_loc=Int(nb_locations(search)/4))
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

![example 1](../../assets/examples/example1.png)

### Insert missing hydrogens

Example: Insert missing H atoms in IRMOF-1

[Input Files](../../../assets/examples/example2.zip)

```jldoctest; output=false
parent = Crystal("IRMOF-1_noH.cif")
infer_bonds!(parent, true)
query = moiety("1,4-C-phenylene_noH.xyz")
replacement = moiety("1,4-C-phenylene.xyz")
repaired_xtal = replace(parent, query => replacement, rand_all=true)
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

![example 2](../../assets/examples/example2.png)

### Repair Disorder and Remove Adsorbates

Example: correct the crystallographic disorder of the PyC-2 ligands and remove guest molecules from the pores.

[Input Files](../../../assets/examples/example3.zip)


```jldoctest; output=false
parent = Crystal("EMEHUB_C2H2.cif", remove_duplicates=true, check_overlap=false)
infer_bonds!(parent, true)
q = moiety("disordered_ligand!.xyz")
p = moiety("4-pyridyl.xyz")
repaired = replace(parent, q => p, rand_all=true)
search = substructure_search(moiety("acetylene.xyz"), repaired, disconnected_component=true)
active = substructure_replace(search, nothing, rand_all=true)
# output
┌ Info: Crystal EMEHUB_C2H2.cif has I 4/m m m space group. I am converting it to P1 symmetry.
└         To prevent this, pass `convert_to_p1=false` to the `Crystal` constructor.
┌ Warning: carbon atom 1 in EMEHUB_C2H2.cif is bonded to more than four atoms!
└ @ Xtals ~/.julia/packages/Xtals/Kf4en/src/bonds.jl:407
┌ Warning: carbon atom 6 in disordered_ligand!.xyz is bonded to more than four atoms!
└ @ Xtals ~/.julia/packages/Xtals/Kf4en/src/bonds.jl:407
┌ Warning: carbon atom 1 in disordered_ligand!.xyz is bonded to more than four atoms!
└ @ Xtals ~/.julia/packages/Xtals/Kf4en/src/bonds.jl:407
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

![example 3](../../assets/examples/example3.png)

### Generate Missing-Linker Defect

Example: create a new channel in UiO-66 via missing-linker defects and formate ion capping.

[Input Files](../../../assets/examples/example4.zip)

```jldoctest; output=false
parent = Crystal("UiO-66.cif")
infer_bonds!(parent, true)
query = moiety("BDC.xyz")
replacement = moiety("formate_caps.xyz")
with_defect = replace(parent, query => replacement; loc=[1, 3, 9])
# output
Name: new_xtal
Bravais unit cell of a crystal.
	Unit cell angles α = 60.000000 deg. β = 60.000000 deg. γ = 60.000000 deg.
	Unit cell dimensions a = 29.634800 Å. b = 29.634800 Å, c = 29.634800 Å
	Volume of unit cell: 18403.100761 Å³

	# atoms = 888
	# charges = 0
	chemical formula: Dict(:Zr => 24, :H => 109, :O => 128, :C => 183)
	space Group: P1
	symmetry Operations:
		'x, y, z'
```

![example 4](../../assets/examples/example4.png)
