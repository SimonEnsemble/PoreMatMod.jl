### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ 968157d0-c8e9-4906-a738-e91f0785c23b
begin
	fig = "mtv-mof"
	figpath = joinpath("/home/adrian/MOFun.jl/figures", fig)
	push!(LOAD_PATH, "/home/adrian/MOFun.jl/src")
	using MOFun, Xtals
	set_path_to_moieties(figpath)
	set_path_to_crystals(figpath)
end

# ╔═╡ a5044e4c-e9e6-4013-917d-3267e9b510b2
begin
	xtal = Crystal("IRMOF-1.cif", remove_duplicates=true)
	strip_numbers_from_atom_labels!(xtal)
	infer_bonds!(xtal, true)
	xtal
end

# ╔═╡ dec571dd-c21b-46d9-bee3-c9d1417946b0
phenylene = moiety("2-!-p-phenylene")

# ╔═╡ b9bb7221-9b10-4026-9f34-f6953f5a03c1
fluorophenylene = moiety("2-F-p-phenylene")

# ╔═╡ 1d4db96c-f638-4d2e-9cfa-1d4fb8f27b07
acetylamidophenylene = moiety("2-acetylamido-p-phenylene")

# ╔═╡ 542057c3-50d4-490f-9adf-349e43c4dba0
intermediate = find_replace(phenylene ∈ xtal, fluorophenylene, nb_loc=8)

# ╔═╡ 7fb483a2-9559-4329-8948-c6bd8a64f5fe
new_xtal = find_replace(phenylene ∈ intermediate, acetylamidophenylene, nb_loc=8)

# ╔═╡ bc1cf01c-6416-4191-b1dc-730fda8419c4
write_cif(new_xtal, "$figpath/new_xtal.cif")

# ╔═╡ Cell order:
# ╠═968157d0-c8e9-4906-a738-e91f0785c23b
# ╠═a5044e4c-e9e6-4013-917d-3267e9b510b2
# ╠═dec571dd-c21b-46d9-bee3-c9d1417946b0
# ╠═b9bb7221-9b10-4026-9f34-f6953f5a03c1
# ╠═1d4db96c-f638-4d2e-9cfa-1d4fb8f27b07
# ╠═542057c3-50d4-490f-9adf-349e43c4dba0
# ╠═7fb483a2-9559-4329-8948-c6bd8a64f5fe
# ╠═bc1cf01c-6416-4191-b1dc-730fda8419c4
