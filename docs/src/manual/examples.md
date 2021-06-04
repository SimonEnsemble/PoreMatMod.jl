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

```
