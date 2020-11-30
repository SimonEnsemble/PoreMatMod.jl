## Examples

### Generate hypothetical structures

To create novel derivatives of a `Crystal`, search it for one of its substructures
and replace with a derivatized moiety.

Example: *ortho* substitution with an acetylamido group at one quarter of the
*p*-phenylene moieties in IRMOF-1.

[[IRMOF-1_clean.cif](../../../assets/IRMOF-1_clean.cif)]
[[2-!-p-phenylene.xyz](../../../assets/2-!-p-phenylene.xyz)]
[[2-acetylamido-p-phenylene.xyz](../../../assets/2-acetylamido-p-phenylene.xyz)]

```julia
xtal = Crystal("IRMOF-1.cif", infer_bonds=:cordero, periodic_boundaries=true)
s_moty = moiety("2-!-p-phenylene")
r_moty = moiety("2-acetylamido-p-phenylene")
search = s_moty ∈ xtal
new_xtal = find_replace(search, nb_loc=nb_locations(search)/4)
```

### Remove solvent molecules

To remove moieties from a `Crystal`, use a `.xyz` file with 0 atoms as the
replacement moiety.

Example: MOF activation.  Remove cyclohexane from an IRMOF-1 unit cell.

[[IRMOF-1_clean.cif](../../../assets/IRMOF-1_clean.cif)]
[[cyclohexane.xyz](../../../assets/cyclohexane.xyz)]
[[nothing.xyz](../../../assets/nothing.xyz)]

```julia
xtal = Crystal("IRMOF-1_clean.cif")
infer_bonds!(xtal, true)
s_moty = moiety("cyclohexane")
r_moty = moiety("nothing")
activated_xtal = find_replace(s_moty ∈ xtal, r_moty, rand_all=true)
```

### Repair disorder

To remove disordered groups of atoms and replace them with ordered moieties,
extract the disordered atoms' coordinates to use as the search moiety, and a
manually corrected or *de novo* structure as the replacement moiety.

Example: Correct the rotational disorder of DABCO linkers in the MOF
ZmID.

[[ZmID.cif](../../../assets/ZmID.cif)]
[[disordered_dabco.xyz](../../../assets/disordered_dabco.xyz)]
[[dabco.xyz](../../../assets/dabco.xyz)]

```julia
# need to set check_overlap due to disordered dabco
xtal = Crystal("ZmID.cif", check_overlap=false)
infer_bonds!(xtal, true)
s_moty = moiety("disordered_dabco")
r_moty = moiety("dabco")
repaired_xtal = find_replace(s_moty ∈ xtal, r_moty, rand_all=true)
```

### Insert missing hydrogens

To correct defects like missing atoms, use the affected substructure as the search
moiety and a manually corrected copy as the replacement moiety.

Example: Insert missing H atoms in IRMOF-1

[[IRMOF-1_noH.cif](../../../assets/IRMOF-1_noH.cif)]
[[p-phenylene_noH.xyz](../../../assets/p-phenylene_noH.xyz)]
[[p-phenylene](../../../assets/p-phenylene.xyz)]

```julia
xtal = Crystal("IRMOF-1_noH.cif")
s_moty = moiety("p-phenylene_noH")
r_moty = moiety("p-phenylene")
repaired_xtal = find_replace(s_moty ∈ xtal, r_moty, rand_all=true)
```
