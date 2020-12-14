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

```julia
xtal = Crystal("IRMOF-1.cif", infer_bonds=:cordero, periodic_boundaries=true)
s_moty = moiety("2-!-p-phenylene")
r_moty = moiety("2-acetylamido-p-phenylene")
search = s_moty ∈ xtal
new_xtal = find_replace(search, nb_loc=nb_locations(search)/4)
```

### Remove solvent molecules

[![example 2](../../../assets/solventexample.png)
Pluto Notebook](../../../examples/solvent.jl)

Remove moieties from a `Crystal`.

Example: MOF activation.  Remove cyclohexane from an IRMOF-1 unit cell.

```julia
xtal = Crystal("IRMOF-1_clean.cif")
infer_bonds!(xtal, true)
s_moty = moiety("cyclohexane")
r_moty = moiety(nothing)
activated_xtal = find_replace(s_moty ∈ xtal, r_moty, rand_all=true)
```

### Repair disorder

[![example 3](../../../assets/disorderexample.png)
Pluto Notebook](../../../examples/disorder.jl)

Remove disordered groups of atoms and replace them with ordered moieties.
Extract the disordered atoms' coordinates to use as the search moiety, and a
manually corrected or *de novo* structure as the replacement moiety.

Example: Correct the rotational disorder of DABCO linkers in the MOF ZmID.

```julia
# need to set check_overlap due to disordered dabco
xtal = Crystal("ZmID.cif", check_overlap=false)
infer_bonds!(xtal, true)
s_moty = moiety("disordered_dabco")
r_moty = moiety("dabco")
repaired_xtal = find_replace(s_moty ∈ xtal, r_moty, rand_all=true)
```

### Insert missing hydrogens

[![example 4](../../../assets/missingHexample.png)
Pluto Notebook](../../../examples/missingH.jl)

To correct defects like missing atoms, use the affected substructure as the search
moiety and a manually corrected copy as the replacement moiety.

Example: Insert missing H atoms in IRMOF-1


```julia
xtal = Crystal("IRMOF-1_noH.cif")
s_moty = moiety("p-phenylene_noH")
r_moty = moiety("p-phenylene")
repaired_xtal = find_replace(s_moty ∈ xtal, r_moty, rand_all=true)
```

### Multiple Transformations

[![example 5](../../../assets/landingpageexample.png)
Pluto Notebook](../../../examples/landingpage.jl)

The example on the [landing page](../../../index.md): repair, activate, and functionalize.

```julia
using MOFun
xtal = Crystal("EMEHUB.cif", infer_bonds=:voronoi, periodic_bonds=true)
repaired = (moiety("disordered_bipy!") => moiety("discrete")) ∈ xtal
active = (moiety("acetylene") => moiety(nothing)) ∈ repaired
novel = (moiety("2-H!-PyC2") => moiety("2-Me-PyC2")) ∈ active
```
