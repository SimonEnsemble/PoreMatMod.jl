### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 5b182574-5656-4035-85c1-89d1c92ff719
begin
	import Pkg
	Pkg.add(url="https://github.com/SimonEnsemble/PoreMatMod.jl")
end

# ╔═╡ 1a73bdc4-9394-4212-9ae8-3b8654a496c0
using PoreMatMod

# ╔═╡ f05cdf05-1661-41c7-863d-15a436791ac4
using PoreMatMod.ExampleHelpers

# ╔═╡ 0b644231-c25a-4d9f-895d-16fc612613ec
md"# Azepine"

# ╔═╡ 916b9bc9-c6c9-444b-9281-e3709c248b75
begin
	benzodiazepine_parent = moiety("benzodiazepine.xyz")
	remove_bonds!(benzodiazepine_parent)
	local new_box = replicate(unit_cube(), (20, 20, 20))
	local new_atoms = Frac(Cart(benzodiazepine_parent.atoms, benzodiazepine_parent.box), new_box)
	local new_charges = Frac(Cart(benzodiazepine_parent.charges, benzodiazepine_parent.box), new_box)
	benzodiazepine_parent = Crystal(benzodiazepine_parent.name, new_box, new_atoms, new_charges)
	infer_bonds!(benzodiazepine_parent, false)
	benzo_fragment = moiety("benzo_fragment.xyz")
	chloro_benzo_fragment = moiety("chloro_benzo_fragment.xyz")
	diazepine_fragment = moiety("diazepine_fragment.xyz")
	N_methyl_diazepine = moiety("N_methyl_diazepine.xyz")
end

# ╔═╡ da1f22fb-8e16-439b-934a-280ed31635df
begin
	intermediate = replace(benzodiazepine_parent, diazepine_fragment => N_methyl_diazepine)
	diazepam = replace(intermediate, benzo_fragment => chloro_benzo_fragment)
end

# ╔═╡ 7d73fa2e-9ee5-4002-97a6-496f4d186512
begin
	local temp = deepcopy(intermediate)
	translate_by!(temp.atoms.coords, Frac([0.5, 0.5, 0.5]))
	wrap!(temp)
	write_cif(temp, "intermediate.cif")
	
	local temp = deepcopy(diazepam)
	translate_by!(temp.atoms.coords, Frac([0.5, 0.5, 0.5]))
	wrap!(temp)
	write_cif(temp, "diazepam.cif")
end

# ╔═╡ 6fe41bf5-cf36-4313-8886-cfca2825c6cd
md"# Slab"

# ╔═╡ faff7bea-0d3a-45fe-8a88-b288b390e35b
function read_poscar(filename::String)::Crystal
	filedata = readlines(filename)
	box_matrix = Matrix(reduce(hcat, [parse.(Float64, row) for row in split.(filedata[3:5])])')
	species = Symbol.(split(filedata[6]))
	atom_counts = parse.(Int, split(filedata[7]))
	nb_atoms = sum(atom_counts)
	species_vec = reduce(vcat, [[element for _ in 1:atom_counts[i]] for (i, element) in enumerate(species)])
	
	box = Box(box_matrix)
	coords = Frac(reduce(hcat, [parse.(Float64, row) for row in split.(filedata[9:nb_atoms+8])]))
	atoms = Atoms(nb_atoms, species_vec, coords)
	charges = Charges(nb_atoms, zeros(nb_atoms), coords)
	
	return Crystal(filename, box, atoms, charges)
end

# ╔═╡ db771098-8097-4748-85c7-ece4faed4e43
begin
	poscar = read_poscar("data/crystals/POSCAR")
	infer_bonds!(poscar, true)
end

# ╔═╡ 35c1c1d5-ea01-4d15-80c3-63330744b037
write_cif(poscar, "poscar.cif")

# ╔═╡ b09ea1c5-67b5-4953-8fef-c5144f09186b
hydrated_Pd2 = moiety("hydrated_Pd2.xyz")

# ╔═╡ b7a6a162-39cc-4137-b803-94c2a5326b6e
OA_Pd2 = moiety("OA_Pd2.xyz")

# ╔═╡ e232cfc8-a954-4089-be59-a9fa5b3a2736
oxidative_addition = replace(poscar, hydrated_Pd2 => OA_Pd2)

# ╔═╡ 908ab709-c12a-455b-994b-89e600d47b06
write_cif(oxidative_addition, "oxidative_addition.cif")

# ╔═╡ Cell order:
# ╟─5b182574-5656-4035-85c1-89d1c92ff719
# ╠═1a73bdc4-9394-4212-9ae8-3b8654a496c0
# ╠═f05cdf05-1661-41c7-863d-15a436791ac4
# ╟─0b644231-c25a-4d9f-895d-16fc612613ec
# ╠═916b9bc9-c6c9-444b-9281-e3709c248b75
# ╠═da1f22fb-8e16-439b-934a-280ed31635df
# ╠═7d73fa2e-9ee5-4002-97a6-496f4d186512
# ╟─6fe41bf5-cf36-4313-8886-cfca2825c6cd
# ╠═faff7bea-0d3a-45fe-8a88-b288b390e35b
# ╠═db771098-8097-4748-85c7-ece4faed4e43
# ╠═35c1c1d5-ea01-4d15-80c3-63330744b037
# ╠═b09ea1c5-67b5-4953-8fef-c5144f09186b
# ╠═b7a6a162-39cc-4137-b803-94c2a5326b6e
# ╠═e232cfc8-a954-4089-be59-a9fa5b3a2736
# ╠═908ab709-c12a-455b-994b-89e600d47b06
