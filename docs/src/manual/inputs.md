```@meta
DocTestSetup = quote
    using PoreMatMod
end
```

# Input Files

This section details the handling of data paths, formatting of input files, and
loading of data into `PoreMatMod`.

## Data Directories

`PoreMatMod` draws its basic working data (atomic masses, covalent radii, etc.) from
[`Xtals.jl`](https://github.com/SimonEnsemble/Xtals.jl/).

Structural data are drawn from `rc[:paths][:crystals]` and `rc[:paths][:moieties]`.
These paths are set at module load time, and default to `./data/crystals` and `./data/moieties`, respectively.

## Input Files and Formats

### Crystals

`PoreMatMod` requires chemical structural data as input.  The first necessary input is a `.cif` or `.cssr` file containing 
atomic coordinates and unit cell information.

The file must be located in `rc[:paths][:crystals]` as described above. In the case of our guiding example, the 
functionalization of IRMOF-1, this means we need to either put 
[`IRMOF-1.cif`](https://raw.githubusercontent.com/SimonEnsemble/PoreMatMod.jl/master/test/data/crystals/IRMOF-1.cif?token=AD3TMGFZCE4WX3J4TDH2BSDAYMO2K) 
into `./data/crystals` or set `rc[:paths][:crystals]` to point `PoreMatMod` to where 
[`IRMOF-1.cif`](https://raw.githubusercontent.com/SimonEnsemble/PoreMatMod.jl/master/test/data/crystals/IRMOF-1.cif?token=AD3TMGFZCE4WX3J4TDH2BSDAYMO2K) 
is located.

### Fragments

The next required input is a `.xyz` file containing the atomic coordinates of a search moiety--a chemical substructure 
to identify in the crystal.  The `.xyz` format is simple: the first line gives the number of input lines which follow,
and each subsequent input line consists of the atom label in the first space-delimited column, followed by 3 columns 
for the atom's Cartesian coordinates in Ångströms.

By default, for use with `PoreMatMod.jl`, `.xyz` data must have clean atom labels, meaning only plain atomic symbols. The 
exception is the use of `!` for indicating atoms which will be altered in a [`replace` operation](../../replace). 
For [substructure searches](../../find) using [`substructure_search`], any `!` tags are ignored (the atoms are 
treated as normal).

The `.xyz` file must be located at `rc[:paths][:moieties]`. For what we want to do with IRMOF-1, the best choice is to 
search for the 
[`*p*-phenylene`](https://raw.githubusercontent.com/SimonEnsemble/PoreMatMod.jl/master/test/data/moieties/p-phenylene.xyz?token=AD3TMGFBEFHHR3NUT4UAP3TAYMPLI) 
moiety that is the core of the BDC linker.

## Loading Files

Load the 
[parent crystal](https://raw.githubusercontent.com/SimonEnsemble/PoreMatMod.jl/master/test/data/crystals/IRMOF-1.cif?token=AD3TMGFZCE4WX3J4TDH2BSDAYMO2K) 
and build the bonding network:

```jldoctest
xtal = Crystal("IRMOF-1.cif")
infer_bonds!(xtal, true)
# output
true
```

[`Crystal`](https://simonensemble.github.io/Xtals.jl/dev/crystal/#Xtals.Crystal) is inherited and re-exported from `Xtals.jl`.
See the [`docs`](https://simonensemble.github.io/Xtals.jl/dev/crystal/#Xtals.Crystal) for more information.

Load the 
[search moiety](https://raw.githubusercontent.com/SimonEnsemble/PoreMatMod.jl/master/test/data/moieties/p-phenylene.xyz?token=AD3TMGFBEFHHR3NUT4UAP3TAYMPLI)
:

```jldoctest
s_moty = moiety("p-phenylene")
# output
Name: p-phenylene
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

Both `xtal` and `s_moty` are `Crystal` objects.

## Documentation

```@docs
moiety
```
