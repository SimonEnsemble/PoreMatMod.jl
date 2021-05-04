# Input Files

This section details the handling of data paths, formatting of input files, and
loading of data into `MOFun`.

## Data Directories

`MOFun` draws its basic working data (atomic masses, covalent radii, etc.) from
[`Xtals.jl`](https://SimonEnsemble.github.io/Xtals.jl/) using `Xtals.PATH_TO_DATA`.
This path is `./data` by default, and can be changed with [`set_path_to_data`](@ref).

```julia
set_path_to_data("path/to/data")
```

Structural data are drawn from `Xtals.PATH_TO_CRYSTALS` and `MOFun.PATH_TO_MOIETIES`.
These paths default to `./data/crystals` and `./data/moieties`, respectively,
and can be changed with [`set_path_to_crystals`](@rwef) and [`set_path_to_moieties`](@ref).

```julia
set_path_to_crystals("path/to/crystals")
set_path_to_moieties("path/to/moieties")
```

The three paths can be set independently of one another, or, if your data are all
in one place, with crystals and moieties in the appropriate sub-directories,
[`set_path_to_data`](@ref) can handle it all.

```julia
# change the data path to "foo" and make the crystals path "foo/crystals"
set_path_to_data("foo", relpath_xtals=true)
# change the data path to "bar" and make the moieties path "bar/crystals"
set_path_to_data("bar", relpath_motys=true)
# change the data paths to "baz/", "baz/crystals", and "baz/moieties"
set_path_to_data("baz", relpaths=true)
```

See the current paths with [`print_file_paths`](@ref).

## Input Files and Formats

### .cif

`MOFun` requires chemical structural data as input.  The first necessary input
is a `.cif` file containing the atomic coordinates and unit cell information for
a crystalline material.

The `.cif` file must be located in `PATH_TO_CRYSTALS` as described above. In the
case of our guiding example, the functionalization of IRMOF-1, this means we need
to either put [`IRMOF-1.cif`](assets/IRMOF-1.cif) into `PATH_TO_DATA/crystals` or
use [`set_path_to_crystals`](@ref) to point `MOFun` to where [`IRMOF-1.cif`](assets/IRMOF-1.cif)
is located.

### .xyz

The next required input is a `.xyz` file containing the atomic coordinates of a
search moiety--a chemical substructure to identify in the crystal.  The `.xyz`
format is simple: the first line gives the number of input lines which follow,
and each subsequent input line consists of the atom label in the first space-
delimited column, followed by 3 columns for the atom's Cartesian coordinates
in Ångströms.

By default, for use with `MOFun.jl`, `.xyz` data must have clean atom labels,
meaning only plain atomic symbols. Exceptions are specifying carbon atom
hybridizations, e.g. `C_sp3`, and manganese/iron/cobalt atom spin states,
e.g. `Mn_lo`; and the use of `!` for indicating atoms which will be altered in
a [`replace` operation](manual/replace). For substructure searches using
[`substructure_search`](@ref), any `!` tags are ignored (the atoms are treated as
normal).

The `.xyz` file must be located in the `moieties` data path. For what we want to
do with IRMOF-1, the best choice is to search for the [`*p*-phenylene`](assets/p-phenylene.xyz)
moiety that is the core of the BDC linker.

## Loading Files

Load the [parent crystal](../../../assets/IRMOF-1.cif):

```julia
xtal = Crystal("IRMOF-1.cif", infer_bonds=:cordero, periodic_boundaries=true)
```

[`Crystal`](https://simonensemble.github.io/Xtals.jl/dev/crystal/#Xtals.Crystal)
is inherited and re-exported from `Xtals.jl`.
See the [`docs`](https://simonensemble.github.io/Xtals.jl/dev/crystal/#Xtals.Crystal)
for more information.

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

```@docs
set_path_to_data
set_path_to_moieties
print_file_paths
moiety
```
