```@meta
DocTestSetup = quote
    using PoreMatMod
end
```

# Reading data from crystal structure and chemical fragment files into `PoreMatMod.jl`

This section details how to load data into `PoreMatMod.jl`, including the handling of paths to data and input file formats.

## Crystal structures

Accepted file formats for crystal structures (containing atomic coordinates and unit cell information) are `.cif` and `.cssr`. See [here](https://en.wikipedia.org/wiki/Crystallographic_Information_File) for information about CIF files.

Crystal structure files (`.cif`, `.cssr`) are read from the path `rc[:paths][:crystals]`.

!!! example
    Read in the crystal structure of [IRMOF-1.cif](../../../assets/inputs/IRMOF-1.cif) and infer its bonding network:

    ```jldoctest; output=false
    parent = Crystal("IRMOF-1.cif")
    infer_bonds!(parent, true) # true b/c we want periodic bonds included
    # output
    true
    ```

the `Crystal` constructor returns a [`Crystal`](@ref) data structure, inherited from `Xtals.jl` (see the [`docs`](https://simonensemble.github.io/Xtals.jl/dev/crystal/#Xtals.Crystal)).

## Query and Replacement Fragments

Accepted file formats for chemical fragments (list of atoms and their Cartesian coordinates) are `.xyz`. See [here](https://en.wikipedia.org/wiki/XYZ_file_format) for information about XYZ files. 

Query and replacement fragment files (`.xyz`) are read from the path `rc[:paths][:moieties]`.

N.b. masked atoms of query fragments must be labeled with `!` for [`replace` operations](../../replace). For [substructure searches](../../find) using `substructure_search`, any `!` tags are ignored (the atoms are treated according to their chemical species).

!!! example
    Read in the chemical fragment [`p-phenylene.xyz`](../../../assets/inputs/p-phenylene.xyz):

    ```jldoctest; output=false
    query = moiety("p-phenylene.xyz")
    # output
    Name: p-phenylene.xyz
    Bravais unit cell of a crystal.
    	Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
    	Unit cell dimensions a = 1.000000 Å. b = 1.000000 Å, c = 1.000000 Å
    	Volume of unit cell: 1.000000 Å³

    	# atoms = 10
    	# charges = 0
    	chemical formula: Dict(:H => 2, :C => 3)
    	space Group: P1
    	symmetry Operations:
    		'x, y, z'
    ```

the [`moiety`](@ref) reader also returns a `Crystal` data structure but with an arbitrary unit cube unit cell.

## Changing the Data Directories

`rc[:paths][:crystals]` and `rc[:paths][:moieties]` default to `./data/crystals` and `./data/moieties`, respectively.

change the paths to where the files are read from by changing `rc[:paths][:crystals]` and `rc[:paths][:moieties]`.

!!! example
    Suppose we wish to store our `.cif` files in `~/my_xtals` and our `.xyz` files in our present working directory.

    ```julia
    rc[:paths][:crystals] = joinpath(homedir(), "my_xtals_folder")
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
