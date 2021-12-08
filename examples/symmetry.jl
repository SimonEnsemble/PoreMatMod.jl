### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ cdf8f316-54a1-11ec-0e6d-b1a485ca5fe8
begin
	import Pkg
	Pkg.develop("PoreMatMod")
end

# ╔═╡ 3bb3966e-a7dd-46c1-b649-33169ce424d2
using PoreMatMod, Bio3DView, PlutoUI, Graphs, LinearAlgebra, StatsBase

# ╔═╡ d3eef6f4-ff15-47f3-8686-2c0cb0fb882d
using PoreMatMod.ExampleHelpers

# ╔═╡ e561a0f6-eaad-4ac8-b814-c85ed3b92a09
using PoreMatMod: Alignment, get_r2p_alignment, aligned_replacement

# ╔═╡ b289f607-6b82-4320-b326-533b982ec6c6
function conglomerate!(parent_substructure::Crystal)
	ref_atom = 1
	# todo
	return
end

# ╔═╡ 536e88c0-3924-4cfb-a4f9-e5ad8c76ea88
begin
	parent = Crystal("NiPyC_fragment_trouble.cif")
	parent = replicate(parent, (2,2,2))
	infer_bonds!(parent, true)
	query = moiety("PyC.xyz")
	replacement = moiety("PyC-CH3.xyz")
	search = query in parent
end

# ╔═╡ 4f8868ec-d566-4525-af7e-c54801e3d7b9
begin
	# arguments
	isomorphisms = search.isomorphisms
	loc_id = 1
	ori_id = 1
end

# ╔═╡ ee7747cc-d5cf-46f3-807b-d3768621d11c
q2p = isomorphisms[loc_id][ori_id]

# ╔═╡ cc9d44ed-2b90-48df-83e9-7181ca2f8286
nb_not_masked = sum(.! occursin.(rc[:r_tag], String.(query.atoms.species)))

# ╔═╡ 0dfe1c41-cadd-4d56-bdf2-106ac17e8df4
@assert sum(occursin.(rc[:r_tag], String.(query.atoms.species))) + nb_not_masked == query.atoms.n

# ╔═╡ 690427dc-c25b-4f23-af11-3d47c132a6b5
@assert ! any(occursin.(rc[:r_tag], String.(query.atoms.species))[1:nb_not_masked]) "all masked atoms should be at the end"

# ╔═╡ 21a25153-9ce8-419f-8c3d-69cdeb0a0d67
q_unmasked_in_r = query[1:nb_not_masked] ∈ replacement

# ╔═╡ a47dc0bd-73a9-49da-9253-9fc28b78be9f
q2r = Dict([q => q_unmasked_in_r.isomorphisms[1][1][q] for q = 1:nb_not_masked])

# ╔═╡ 96f583de-ff3b-4287-bf16-132eaf9adc05
r2p = Dict([r => q2p[q] for (q, r) in q2r])

# ╔═╡ 04e61f68-7336-4c1d-888b-50fba2cae00d
r2p_alignment = get_r2p_alignment(replacement, parent, r2p)

# ╔═╡ 0cbf3baa-24fc-4022-a08e-8b24379d1500
replacement_to_install = aligned_replacement(replacement, parent, r2p_alignment)

# ╔═╡ 2b342469-8967-4817-9439-ea6625321c41
trim(parent::Crystal, q2p::Dict{Int, Int}) = parent[[p for p = 1:parent.atoms.n if ! (p in values(q2p))]]

# ╔═╡ 586f38c8-4bd5-4cd5-a748-896f795afc6f
child = trim(parent, q2p) + replacement_to_install

# ╔═╡ d8ea8d5e-f17a-412b-8461-15ba6d9621ec
write_cif(child)

# ╔═╡ eae2225c-40f0-4d68-a9a2-43a39a82f029
view_structure(child)

# ╔═╡ Cell order:
# ╠═cdf8f316-54a1-11ec-0e6d-b1a485ca5fe8
# ╠═3bb3966e-a7dd-46c1-b649-33169ce424d2
# ╠═d3eef6f4-ff15-47f3-8686-2c0cb0fb882d
# ╠═e561a0f6-eaad-4ac8-b814-c85ed3b92a09
# ╠═b289f607-6b82-4320-b326-533b982ec6c6
# ╠═536e88c0-3924-4cfb-a4f9-e5ad8c76ea88
# ╠═4f8868ec-d566-4525-af7e-c54801e3d7b9
# ╠═ee7747cc-d5cf-46f3-807b-d3768621d11c
# ╠═cc9d44ed-2b90-48df-83e9-7181ca2f8286
# ╠═0dfe1c41-cadd-4d56-bdf2-106ac17e8df4
# ╠═690427dc-c25b-4f23-af11-3d47c132a6b5
# ╠═21a25153-9ce8-419f-8c3d-69cdeb0a0d67
# ╠═a47dc0bd-73a9-49da-9253-9fc28b78be9f
# ╠═96f583de-ff3b-4287-bf16-132eaf9adc05
# ╠═04e61f68-7336-4c1d-888b-50fba2cae00d
# ╠═0cbf3baa-24fc-4022-a08e-8b24379d1500
# ╠═2b342469-8967-4817-9439-ea6625321c41
# ╠═586f38c8-4bd5-4cd5-a748-896f795afc6f
# ╠═d8ea8d5e-f17a-412b-8461-15ba6d9621ec
# ╠═eae2225c-40f0-4d68-a9a2-43a39a82f029
