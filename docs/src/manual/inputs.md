# Preparing Inputs

## Input Files

`MOFun.jl` requires chemical structural data as input.  The first necessary input
is a `.cif` file containing the atomic coordinates and unit cell information for
a crystalline material.

The `.cif` file must be located in the `crystals` subdirectory of
the data path for `Xtals.jl`.  Set the data paths with:

```julia
set_path_to_data("path to data", relpath_xtals=true, relpath_motys=true)
```

See the [Xtals.jl docs](https://github.com/SimonEnsemble/Xtals.jl) for more
information on setting the paths.

The next required input is a `.xyz` file containing the atomic coordinates of a
search moiety--a chemical substructure to identify in the crystal.

The `.xyz`
file must be located in the `moieties` data path for `MOFun.jl`.  Set the data
path with:

```julia
set_path_to_moieties("path to moieties")
```

## Loading Files

Load the [parent crystal](../../../assets/IRMOF-1.cif):

```julia
xtal = Crystal("IRMOF-1.cif")
```

The crystal needs to have clean site labels for bond inference.  Clean up the
labels as necessary and infer the bonds:

```julia
strip_numbers_from_atom_labels!(xtal) # if needed
infer_bonds!(xtal, true) # true -> bond over periodic boundaries
```

Load the [search moiety](../../../assets/p-phenylene.xyz):

```julia
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
[Xtals docs](https://simonensemble.github.io/Xtals.jl/dev/crystal/#Xtals.Crystal)

```@docs
set_path_to_data
set_path_to_moieties
print_file_paths
moiety
```
