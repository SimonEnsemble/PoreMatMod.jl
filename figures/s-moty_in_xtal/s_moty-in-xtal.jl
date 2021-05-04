### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ 968157d0-c8e9-4906-a738-e91f0785c23b
begin
	push!(LOAD_PATH, "/home/adrian/MOFun.jl/src")
	using MOFun, Xtals
	set_path_to_moieties("/home/adrian/MOFun.jl/figures/s-moty_in_xtal")
	set_path_to_crystals("/home/adrian/MOFun.jl/figures/s-moty_in_xtal")
end

# ╔═╡ a5044e4c-e9e6-4013-917d-3267e9b510b2
begin
	xtal = Crystal("IRMOF-1.cif")
	strip_numbers_from_atom_labels!(xtal)
	infer_bonds!(xtal, true)
	s_moty = moiety("p-phenylene")
	search = s_moty ∈ xtal
end

# ╔═╡ 33ec822d-db8a-4fc9-88b9-eb0a3eb34bb7
search.results

# ╔═╡ fe23ad61-948f-4d72-9310-5f789b60a8be
phenylene_ids = collect(Iterators.flatten([search.results[i].isomorphism[1] for i ∈ 1:search.results.ngroups]))

# ╔═╡ 81ea5d21-12fa-462c-9bed-d8778bcd52bd
phenylenes = xtal[phenylene_ids]

# ╔═╡ cb885730-5fea-49a7-b013-85b5be6d4805
Xtals.write_xyz(phenylenes, joinpath(pwd(), "phenylenes.xyz"))

# ╔═╡ 2652dde8-b1df-484d-b7ae-964ed02dcc1b
Xtals.write_vtk(phenylenes.box, joinpath(pwd(), "phenylenes.vtk"))

# ╔═╡ b657e688-b8b6-43ea-a2cf-00212486890f
write_cif(phenylenes, joinpath(pwd(), "phenylenes.cif"))

# ╔═╡ Cell order:
# ╠═968157d0-c8e9-4906-a738-e91f0785c23b
# ╠═a5044e4c-e9e6-4013-917d-3267e9b510b2
# ╠═33ec822d-db8a-4fc9-88b9-eb0a3eb34bb7
# ╠═fe23ad61-948f-4d72-9310-5f789b60a8be
# ╠═81ea5d21-12fa-462c-9bed-d8778bcd52bd
# ╠═cb885730-5fea-49a7-b013-85b5be6d4805
# ╠═2652dde8-b1df-484d-b7ae-964ed02dcc1b
# ╠═b657e688-b8b6-43ea-a2cf-00212486890f
