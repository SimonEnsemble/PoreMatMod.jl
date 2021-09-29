```@meta
DocTestSetup = quote
    using PoreMatMod
end
```

## Examples

The [Pluto notebooks](https://github.com/fonsp/Pluto.jl) (`.jl`) containing the Julia code and the input files for all examples are in the [examples directory](https://github.com/SimonEnsemble/PoreMatMod.jl/tree/master/examples).

Click on the images or descriptions below to see Pluto notebooks demonstrating each use case.

!!! example "Example: Generate hypothetical structures"
    Decorate IRMOF-1 with a functional group: *ortho* substitution with an acetylamido group at one quarter of the *p*-phenylene moieties.

    [link](../examples/make_hypothetical_MOF.jl.html)
    
    [
        ![example 1](../assets/examples/example1.png)
    ](../../examples/make_hypothetical_MOF.jl.html)


!!! example "Example: Using different replacement modes"
    Replace BDC with nitro-BDC in IRMOF-1 using the different replacement modes to obtain different results.

    [link](../examples/replacement_modes.jl.html)

    [
        ![example1.5](../assets/examples/example1.5.png)
    ](../../examples/replacement_modes.jl.html)

!!! example "Example: Insert missing hydrogen atoms"
     Insert missing H atoms in IRMOF-1.

     [link](../examples/correct_missing_Hs.jl.html)

    [
        ![example 2](../assets/examples/example2.png)
    ](../../examples/correct_missing_Hs.jl.html)


!!! example "Example: Repair Disorder and Remove Adsorbates"
    Correct the crystallographic disorder of the PyC-2 ligands and remove guest molecules from the pores of SIFSIX-Cu-2-i.

    [link](../examples/disorder_and_guests.jl.html)

    [
        ![example 3](../assets/examples/example3.png)
    ](../../examples/disorder_and_guests.jl.html)

!!! example "Example: Generate Missing-Linker Defects"
    Create a new channel in UiO-66 by introducing missing-linker defects and formate ion capping.

    [link](../examples/missing_linker_defect.jl.html)

    [
        ![example 4](../assets/examples/example4.png)
    ](../../examples/missing_linker_defect.jl.html)
