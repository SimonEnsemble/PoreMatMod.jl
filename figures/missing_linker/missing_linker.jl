### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 968157d0-c8e9-4906-a738-e91f0785c23b
begin
	fig = "missing_linker"
	figpath = joinpath("/home/adrian/MOFun.jl/figures", fig)
	push!(LOAD_PATH, "/home/adrian/MOFun.jl/src")
	using MOFun, Xtals
	set_path_to_moieties(figpath)
	set_path_to_crystals(figpath)
end

# ╔═╡ a5044e4c-e9e6-4013-917d-3267e9b510b2
begin
	xtal = replicate(Crystal("RUBTAK01_SL.cif", remove_duplicates=true), (2,2,2))
	strip_numbers_from_atom_labels!(xtal)
	infer_bonds!(xtal, true)
	xtal
end

# ╔═╡ f8c4b802-bbbf-4c9c-b2ce-d4635c869fac
write_cif(xtal, "$figpath/original.cif")

# ╔═╡ f01758a1-f57e-4d8a-8589-7db914bd743c
bdc = moiety("bdc")

# ╔═╡ 4190bb8d-1687-4274-953e-8f474f4025d5
search = substructure_search(bdc, xtal)

# ╔═╡ 542057c3-50d4-490f-9adf-349e43c4dba0
with_defect = find_replace(search, moiety(nothing), loc=[1,3])

# ╔═╡ bc1cf01c-6416-4191-b1dc-730fda8419c4
write_cif(with_defect, "$figpath/with_defect.cif")

# ╔═╡ Cell order:
# ╠═968157d0-c8e9-4906-a738-e91f0785c23b
# ╠═a5044e4c-e9e6-4013-917d-3267e9b510b2
# ╠═f8c4b802-bbbf-4c9c-b2ce-d4635c869fac
# ╠═f01758a1-f57e-4d8a-8589-7db914bd743c
# ╠═4190bb8d-1687-4274-953e-8f474f4025d5
# ╠═542057c3-50d4-490f-9adf-349e43c4dba0
# ╠═bc1cf01c-6416-4191-b1dc-730fda8419c4
