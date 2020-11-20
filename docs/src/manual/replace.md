# Find/Replace Operations

## Inputs

To guide replacement, the `.xyz` file input for the search moiety must be
altered.  Simply adding `!` after the atomic symbol tags that atom for replacement.

In the example of finding and replacing an *ortho* H of the *p*-phenylene moieties of
IRMOF-1, the input should have one H atom tagged, like so:

```
H!         1.06706        0.70670        1.48683
```

Load it in as before.

```
s_moty = moiety("2-!-p-phenylene")
```

A final additional file input is required: a `.xyz` file containing the replacement
moiety.  This structure must contain a subgraph isomorphic to the non-tagged atoms
of the search moiety.

To replace *p*-phenylene moieties with 2-acetylamido-*p*-phenylene moieties, provide
the appropriate data and load the replacement moiety:

```
r_moty = moiety("2-acetylamido-p-phenylene")
```

## Replacement Modes

With all three file inputs loaded (parent crystal `xtal`, search moiety `s_moty`
and replacement moiety `r_moty`) and a substructure search `search` performed,
replacements can be made.

`MOFun.jl` has several replacement modes, one of which must be specified.

### Random orientation at each location

Random configurations will be chosen for each location in `search.results`, so
that each occurrence of the search moiety in the parent crystal is replaced.

```
find_replace(search, r_moty, rand_all=true)
```

### Random orientation at n random locations

The parent crystal's search moieties will be replaced using random configurations
at each of $n$ locations.

```
find_replace(search, r_moty, nb_loc=4)
```

### Random orientation at specific locations

The parent crystal is modified using random configurations at a list of specified
locations.

```
find_replace(search, r_moty, loc=[13, 20])
```

### Specific orientations at specific locations

Specific replacements are made.

```
find_replace(search, r_moty, loc=[13, 20], ori=[1, 1])
```

## Documentation

```@meta
CurrentModule = MOFun
DocTestSetup = quote
    using MOFun
end
```

```@docs
find_replace
```

## Examples

### Generate hypothetical structures

To create novel derivatives of a `Crystal`, search it for one of its substructures
and replace with a derivatized moiety.

Example: *ortho* substitution with an acetylamido group at one quarter of the
*p*-phenylene moieties in IRMOF-1.

[[IRMOF-1_clean.cif](../../../assets/IRMOF-1_clean.cif)]
[[2-!-p-phenylene.xyz](../../../assets/2-!-p-phenylene.xyz)]
[[2-acetylamido-p-phenylene.xyz](../../../assets/2-acetylamido-p-phenylene.xyz)]

```
xtal = Crystal("IRMOF-1_clean.cif")
infer_bonds!(xtal, true)
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

```
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

```
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

```
xtal = Crystal("IRMOF-1_noH.cif")
s_moty = moiety("p-phenylene_noH")
r_moty = moiety("p-phenylene")
repaired_xtal = find_replace(s_moty ∈ xtal, r_moty, rand_all=true)
```
