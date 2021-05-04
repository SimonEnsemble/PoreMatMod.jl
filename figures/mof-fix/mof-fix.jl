### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ 968157d0-c8e9-4906-a738-e91f0785c23b
begin
	fig = "mof-fix"
	figpath = joinpath("/home/adrian/MOFun.jl/figures", fig)
	push!(LOAD_PATH, "/home/adrian/MOFun.jl/src")
	using MOFun, Xtals
	set_path_to_moieties(figpath)
	set_path_to_crystals(figpath)
end

# ╔═╡ a5044e4c-e9e6-4013-917d-3267e9b510b2
begin
	xtal = Crystal("EMEHUB_C2H2.cif", remove_duplicates=true)
	strip_numbers_from_atom_labels!(xtal)
	infer_bonds!(xtal, true)
	xtal
end

# ╔═╡ dec571dd-c21b-46d9-bee3-c9d1417946b0
disordered_ligand = moiety("disordered_ligand!")

# ╔═╡ b9bb7221-9b10-4026-9f34-f6953f5a03c1
pyridyl = moiety("4-pyridyl")

# ╔═╡ dcb49f4a-acf2-4c75-9ed5-28f20ca1b844
no_disorder = (disordered_ligand => pyridyl) ∈ xtal

# ╔═╡ 1f9e5f28-b427-44ff-9aca-16672865becf
write_cif(no_disorder, "$figpath/no_disorder.cif")

# ╔═╡ f01758a1-f57e-4d8a-8589-7db914bd743c
acetylene = moiety("C2H2!")

# ╔═╡ 4190bb8d-1687-4274-953e-8f474f4025d5
search = substructure_search(acetylene, no_disorder, exact=true)

# ╔═╡ 542057c3-50d4-490f-9adf-349e43c4dba0
clean_xtal = find_replace(search, moiety(nothing), rand_all=true)

# ╔═╡ a9e4efcd-6786-45f6-a7b1-b53601a13f9a
Xtals.write_xyz(clean_xtal, "$figpath/clean_xtal.cif")

# ╔═╡ 3f37a6ee-bfed-4601-b640-c023e4616c19
write_bond_information(clean_xtal, "$figpath/clean_xtal.vtk")

# ╔═╡ bc1cf01c-6416-4191-b1dc-730fda8419c4
write_cif(clean_xtal, "$figpath/clean_xtal.cif")

# ╔═╡ Cell order:
# ╠═968157d0-c8e9-4906-a738-e91f0785c23b
# ╠═a5044e4c-e9e6-4013-917d-3267e9b510b2
# ╠═dec571dd-c21b-46d9-bee3-c9d1417946b0
# ╠═b9bb7221-9b10-4026-9f34-f6953f5a03c1
# ╠═dcb49f4a-acf2-4c75-9ed5-28f20ca1b844
# ╠═1f9e5f28-b427-44ff-9aca-16672865becf
# ╠═f01758a1-f57e-4d8a-8589-7db914bd743c
# ╠═4190bb8d-1687-4274-953e-8f474f4025d5
# ╠═542057c3-50d4-490f-9adf-349e43c4dba0
# ╠═a9e4efcd-6786-45f6-a7b1-b53601a13f9a
# ╠═3f37a6ee-bfed-4601-b640-c023e4616c19
# ╠═bc1cf01c-6416-4191-b1dc-730fda8419c4
