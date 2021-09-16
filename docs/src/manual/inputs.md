```@meta
DocTestSetup = quote
    using PoreMatMod
end
```

# Reading data from crystal structure and chemical fragment files into `PoreMatMod.jl`

This section details how to load data into `PoreMatMod.jl`, including the handling of paths to data and input file formats.

## Crystal structures

Accepted file formats for crystal structures (containing atomic coordinates and unit cell information) are `.cif` and `.cssr`.

Crystal structure files (`.cif`, `.cssr`) are read from the path `rc[:paths][:crystals]`.

_Example_: Read [IRMOF-1.cif](../../../assets/inputs/IRMOF-1.cif) and infer the bonding network:

```jldoctest; output=false
parent_xtal = Crystal("IRMOF-1.cif")
infer_bonds!(parent_xtal, true) # true b/c we want periodic bonds included
# output
true
```

the `Crystal` data structure (`parent_xtal::Crystal`) is inherited from `Xtals.jl` (see the [`docs`](https://simonensemble.github.io/Xtals.jl/dev/crystal/#Xtals.Crystal)).

## Query and Replacement Fragments

Accepted file formats for chemical fragments (list of atoms and their Cartesian coordinates) are `.xyz`.

Query and replacement fragment files (`.xyz`) are read from the path `rc[:paths][:moieties]`.

_Example_: Read [`p-phenylene.xyz`](../../../assets/inputs/p-phenylene.xyz):

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

the `query` is also stored as a `Crystal` with an (arbitrary) unit cube unit cell.

## Changing the Data Directories

`rc[:paths][:crystals]` and `rc[:paths][:moieties]` default to `./data/crystals` and `./data/moieties`, respectively.

change the paths to where the files are read from by changing `rc[:paths][:crystals]` and `rc[:paths][:moieties]`. e.g.:
```julia
rc[:paths][:crystals] = joinpath(homedir(), "my_xtals_folder") # look in ~/my_xtals_folder
rc[:paths][:moiety] = pwd()                                    # look in present working directory
```

## Other data
`PoreMatMod.jl` draws atomic masses and covalent radii from [`Xtals.jl`](https://github.com/SimonEnsemble/Xtals.jl/).


For use with `PoreMatMod.jl`, `.xyz` data must have clean atom labels, meaning only plain atomic symbols. 
The exception is the use of `!` for indicating atoms which will be altered in a [`replace` operation](../../replace). 
For [substructure searches](../../find) using [`substructure_search`], any `!` tags are ignored (the atoms are treated as normal).

![query fragment](../../assets/inputs/query.png)

The `.xyz` file must be located at `rc[:paths][:moieties]`.
For what we want to do with IRMOF-1, the best choice is to search for the [`p-phenylene.xyz`](../../../assets/inputs/p-phenylene.xyz) moiety that is the core of the BDC linker.

## Loading Files



Both `parent` and `query` are `Crystal` objects.

## Documentation

```@docs
moiety
```
