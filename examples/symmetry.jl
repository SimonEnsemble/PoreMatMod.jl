### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ cdf8f316-54a1-11ec-0e6d-b1a485ca5fe8
begin
	import Pkg
	Pkg.develop("PoreMatMod")
	Pkg.develop("Xtals")
end

# ╔═╡ 3bb3966e-a7dd-46c1-b649-33169ce424d2
using PoreMatMod, Bio3DView, PlutoUI, Graphs, LinearAlgebra, StatsBase, MetaGraphs

# ╔═╡ d3eef6f4-ff15-47f3-8686-2c0cb0fb882d
using PoreMatMod.ExampleHelpers

# ╔═╡ e561a0f6-eaad-4ac8-b814-c85ed3b92a09
using PoreMatMod: Alignment, Installation, optimal_replacement

# ╔═╡ 05ab6500-272c-43a9-8ee7-768145ef0c5c
md"## dev stuff"

# ╔═╡ 8e4fd703-f53e-4056-9466-66f07bacad8d
md"## inputs"

# ╔═╡ 536e88c0-3924-4cfb-a4f9-e5ad8c76ea88
begin
	parent = Crystal("NiPyC_fragment_trouble.cif")
	#parent = replicate(parent, (2,2,2))
	infer_bonds!(parent, true)
	query = moiety("PyC.xyz")
	replacement = moiety("PyC-CH3.xyz")
	search = query in parent
end

# ╔═╡ c85aa32e-9c56-4b04-85cc-f4db3174bd77
md"## setup"

# ╔═╡ cc9d44ed-2b90-48df-83e9-7181ca2f8286
begin
	nb_not_masked = sum(.! occursin.(rc[:r_tag], String.(query.atoms.species)))
	q_unmasked_in_r = query[1:nb_not_masked] ∈ replacement
	q2r = Dict([q => q_unmasked_in_r.isomorphisms[1][1][q] for q = 1:nb_not_masked])
	installation = optimal_replacement(search, replacement, q2r, 1)
end

# ╔═╡ 37a4feae-b957-44a5-8631-3811eecfaef5
md"## add crystals together"

# ╔═╡ 94cfc1d5-5ee9-43d5-b936-0541d260f511
function install_replacements(parent::Crystal,
		replacements::Vector{Installation})::Crystal
	child = deepcopy(parent)
	
	obsolete_atoms = Int[] # to delete at the end
	
	# loop over replacements to install
	for installation in replacements
		replacement, q2p, r2p = 
			installation.aligned_replacement, installation.q2p, installation.r2p
		#add into parent
		child = +(child, replacement, check_overlap=false)
		
		# reconstruct bonds
		for (r, p) in r2p # p is in parent_subst
			p_nbrs = neighbors(parent.bonds, p)
			for p_nbr in p_nbrs
				if ! (p_nbr in values(q2p)) # p_nbr not in parent_subst
					# need bond nbr => r in child, where r is in replacement
					e = (p_nbr, child.atoms.n - replacement.atoms.n + r)
					add_edge!(child.bonds, e)
					## TODO make this always true (currently an assumption):
					set_prop!(child.bonds, e[1], e[2], :cross_boundary, 
						get_prop(parent.bonds, p, p_nbr, :cross_boundary)
					)
				end
			end
		end
		
		# accumulate atoms to delete
		obsolete_atoms = vcat(obsolete_atoms, values(q2p)...)
	end
	
	# delete obsolete atoms
	obsolete_atoms = unique(obsolete_atoms)
	keep_atoms = [p for p = 1:child.atoms.n if ! (p in obsolete_atoms)]
	child = child[keep_atoms]
	# return result
	return child
end

# ╔═╡ 586f38c8-4bd5-4cd5-a748-896f795afc6f
child = install_replacements(parent, [installation])

# ╔═╡ f69bf2e7-87f8-403a-ae23-cf4b7bf7d9f1
md"## output"

# ╔═╡ 75e076ad-4cc0-42c3-9f2e-8138a810ae98
@assert child.atoms.n - parent.atoms.n - replacement.atoms.n + query.atoms.n == 0

# ╔═╡ 9976fb79-210d-4eac-a99a-b81e02586039
@assert ne(child.bonds) - ne(parent.bonds) + ne(replacement.bonds) - ne(query.bonds) == 0 ne.([child.bonds, parent.bonds, replacement.bonds, query.bonds])

# ╔═╡ d8ea8d5e-f17a-412b-8461-15ba6d9621ec
write_cif(child)

# ╔═╡ eae2225c-40f0-4d68-a9a2-43a39a82f029
view_structure(child)

# ╔═╡ Cell order:
# ╟─05ab6500-272c-43a9-8ee7-768145ef0c5c
# ╠═cdf8f316-54a1-11ec-0e6d-b1a485ca5fe8
# ╠═3bb3966e-a7dd-46c1-b649-33169ce424d2
# ╠═d3eef6f4-ff15-47f3-8686-2c0cb0fb882d
# ╠═e561a0f6-eaad-4ac8-b814-c85ed3b92a09
# ╟─8e4fd703-f53e-4056-9466-66f07bacad8d
# ╠═536e88c0-3924-4cfb-a4f9-e5ad8c76ea88
# ╟─c85aa32e-9c56-4b04-85cc-f4db3174bd77
# ╠═cc9d44ed-2b90-48df-83e9-7181ca2f8286
# ╟─37a4feae-b957-44a5-8631-3811eecfaef5
# ╠═94cfc1d5-5ee9-43d5-b936-0541d260f511
# ╠═586f38c8-4bd5-4cd5-a748-896f795afc6f
# ╟─f69bf2e7-87f8-403a-ae23-cf4b7bf7d9f1
# ╠═75e076ad-4cc0-42c3-9f2e-8138a810ae98
# ╠═9976fb79-210d-4eac-a99a-b81e02586039
# ╠═d8ea8d5e-f17a-412b-8461-15ba6d9621ec
# ╠═eae2225c-40f0-4d68-a9a2-43a39a82f029
