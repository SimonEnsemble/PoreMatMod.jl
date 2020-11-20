# Preparing Inputs

## Input Files

`MOFun.jl` requires chemical structural data as input.  The first necessary input
is a `.cif` file containing the atomic coordinates and unit cell information for
a crystalline material.

The `.cif` file must be located in the `crystals` subdirectory of
the data path for `PorousMaterials.jl`.  Set the data path with:

```
@eval PorousMaterials PATH_TO_DATA="path to data"
```

The next required input is a `.xyz` file containing the atomic coordinates of a
search moiety--a chemical substructure to identify in the crystal.

The `.xyz`
file must be located in the `moieties` data path for `MOFun.jl`.  Set the data
path with:

```
@eval MOFun PATH_TO_MOIETIES="path to moieties"
```

## Loading Files

Load the [parent crystal](../../../assets/IRMOF-1.cif):

```
xtal = Crystal("IRMOF-1.cif")
```

The crystal needs to have clean site labels for bond inference.  Clean up the
labels as necessary and infer the bonds:

```
strip_numbers_from_atom_labels!(xtal) # if needed
infer_bonds!(xtal, true) # true -> bond over periodic boundaries
```

And load the [search moiety](../../../assets/p-phenylene.xyz):

```
s_moty = moiety("p-phenylene")
```

## Documentation

```@meta
CurrentModule = MOFun
DocTestSetup = quote
    using MOFun
end
```

See the docstring for `Crystal` at the
[PorousMaterials docs](https://simonensemble.github.io/PorousMaterials.jl/dev/crystal/#PorousMaterials.Crystal)

```@docs
moiety
```
