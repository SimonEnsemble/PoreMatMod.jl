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

Note that the order of atoms as stored in the `Crystal` returned by `moiety` may be different than the order in the file.
This is to speed up structure searches.
If it is important that your `Crystal` have its atoms indexed identically to the source file, one solution is to save a new version of the file using [`write_xyz`](@ref).

!!! example

    Sort the atoms in [`glycine_res.xyz`](../../../assets/inputs/glycine_res.xyz):

    ```jldoctest; output=false
    # read the original data, pre-sorting the atoms
    q1 = moiety("glycine_res.xyz")

    # q1 is now indexed differently than the input data

    # save a new source file
    write_xyz(q1, joinpath(rc[:paths][:moieties], "glycine_res_sorted.xyz"))

    # q2 is ordered the same as the new file
    q2 = moiety("glycine_res_sorted.xyz")

    # q1 and q2 are identical
    @assert isapprox(q1.atoms.coords.xf, q2.atoms.coords.xf; atol=0.01)

    # output
    
    ```

The pre-sorting can also be disabled for non-!-tagged atoms, but at the risk of degraded search efficiency.

!!! example

    Load [`glycine_res.xyz`](../../../assets/inputs/glycine_res.xyz) without changing the atom order:

    ```jldoctest; output=false
    moiety("glycine_res.xyz"; presort=false)

    # output
    Crystal(C₂H₃NO, periodic = TTT):
        bounding_box      : [       1        0        0;
                             6.12323e-17        1        0;
                             6.12323e-17 6.12323e-17        1]u"Å"

        Atoms{Frac}(1, [:N], Frac([-2.4029152; -2.23405082; 0.0;;]))
        Atoms{Frac}(1, [:H], Frac([-1.4033551999999998; -2.26371682; 0.0;;]))
        Atoms{Frac}(1, [:C], Frac([-3.0898692; -0.95823882; 0.0;;]))
        Atoms{Frac}(1, [:C], Frac([-2.0853462; 0.18518518; 0.0;;]))
        Atoms{Frac}(1, [:H], Frac([-3.7147022; -0.88143582; -0.889823;;]))
        Atoms{Frac}(1, [:H], Frac([-3.7147022; -0.88143582; 0.889823;;]))
        Atoms{Frac}(1, [:O], Frac([-0.8513402; -0.06139382; 0.0;;]))
    ```

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
