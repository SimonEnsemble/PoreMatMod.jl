```@meta
DocTestSetup = quote
    using PoreMatMod
end
```

## Examples

The code and data files for all examples can be found here: [Example Notebooks and Data](../../../assets/examples/examples.zip)

Click on the images or descriptions below to see Pluto notebooks demonstrating each use case.

## Generate hypothetical structures

[
    Example: *ortho* substitution with an acetylamido group at one quarter of the *p*-phenylene moieties in IRMOF-1.
    ![example 1](../../assets/examples/example1.png)
](../../../examples/make_hypothetical_MOF.jl.html)

## Insert missing hydrogens

[
    Example: Insert missing H atoms in IRMOF-1
    ![example 2](../../assets/examples/example2.png)
](../../../examples/correct_missing_Hs.jl.html)

## Repair Disorder and Remove Adsorbates

[
    Example: correct the crystallographic disorder of the PyC-2 ligands and remove guest molecules from the pores.
    ![example 3](../../assets/examples/example3.png)
](../../../examples/disorder_and_guests.jl.html)

## Generate Missing-Linker Defect

[
    Example: create a new channel in UiO-66 via missing-linker defects and formate ion capping.
    ![example 4](../../assets/examples/example4.png)
](../../../examples/missing_linker_defect.jl.html)
