```@meta
DocTestSetup = quote
    using PoreMatMod
end
```

# Reading data from crystal structure and chemical fragment files into `PoreMatMod.jl`

This section details how to load data into `PoreMatMod.jl`, including the handling of paths to data and input file formats.

## Crystal structures

Accepted file formats for crystal structures (containing atomic coordinates and unit cell information) are `.cif` (see [here](https://en.wikipedia.org/wiki/Crystallographic_Information_File)) and `.cssr`.

Crystal structure files (`.cif`, `.cssr`) are read from the path `rc[:paths][:crystals]`.

!!! example
    
    Read in the crystal structure of [IRMOF-1.cif](../../../assets/inputs/IRMOF-1.cif) and infer its bonding graph:
    
    ```jldoctest; output=false
    parent = Crystal("IRMOF-1.cif")
    infer_bonds!(parent, true) # true b/c we want periodic bonds included
    
    # output
    
    true
    ```

The `Crystal` constructor returns a [`Crystal`](@ref) data structure.
The [`infer_bonds!`](@ref) function infers the bonding graph of the crystal structure (nodes: atoms, edges: bonds) based on interatomic distances---necessary for subgraph matching.
Both `Crystal` and `infer_bonds!` are inherited from `Xtals.jl` (see the [`docs`](https://simonensemble.github.io/Xtals.jl/dev/crystal/#Xtals.Crystal)).

## Query and Replacement Fragments

Accepted file formats for chemical fragments (list of atoms and their Cartesian coordinates) are `.xyz` (see [here](https://en.wikipedia.org/wiki/XYZ_file_format)).

Query and replacement fragment files (`.xyz`) are read from the path `rc[:paths][:moieties]`.

N.b. masked atoms of query fragments must be labeled with `!` for [`replace` operations](../../replace). For [substructure searches](../../find) using `substructure_search`, any `!` tags are ignored (the atoms are treated according to their chemical species).

!!! example
    
    Read in the chemical fragment [`p-phenylene.xyz`](../../../assets/inputs/p-phenylene.xyz):
    
    ```jldoctest; output=false
    query = moiety("p-phenylene.xyz")

    # output

    Crystal(C₆H₄, periodic = TTT):
        bounding_box      : [       1        0        0;
                             6.12323e-17        1        0;
                             6.12323e-17 6.12323e-17        1]u"Å"
    ```

The [`moiety`](@ref) reader also returns a `Crystal` data structure but with a (arbitrary) unit cube unit cell.

## Changing the Data Directories

`rc[:paths][:crystals]` and `rc[:paths][:moieties]` default to `./data/crystals` and `./data/moieties`, respectively.

To change the paths from where the input files are read, change `rc[:paths][:crystals]` and `rc[:paths][:moieties]`.

!!! example
    
    Suppose we wish to store our `.cif` files in `~/my_xtals` and our `.xyz` files in our present working directory.
    
    ```julia
    rc[:paths][:crystals] = joinpath(homedir(), "my_xtals")
    rc[:paths][:moiety] = pwd()
    ```

## Other data

`PoreMatMod.jl` draws atomic masses and covalent radii from [`Xtals.jl`](https://github.com/SimonEnsemble/Xtals.jl/).

## Detailed documentation for functions

```@docs
moiety
Crystal
infer_bonds!
BondingRule
strip_numbers_from_atom_labels!
```
