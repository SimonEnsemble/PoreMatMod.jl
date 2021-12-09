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
using PoreMatMod, Bio3DView, PlutoUI, Graphs, LinearAlgebra, StatsBase, MetaGraphs

# ╔═╡ d3eef6f4-ff15-47f3-8686-2c0cb0fb882d
using PoreMatMod.ExampleHelpers

# ╔═╡ e561a0f6-eaad-4ac8-b814-c85ed3b92a09
using PoreMatMod: Alignment, aligned_replacement

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

# ╔═╡ 67692b25-bc22-4663-be60-8f3a5ee4a4cd
md"## helpers (move to PMM)"

# ╔═╡ b289f607-6b82-4320-b326-533b982ec6c6
function conglomerate!(parent_substructure::Crystal)
	@assert length(connected_components(parent_substructure.bonds)) == 1 "multiple connected components in parent"
	# snip cross-PB bonds to generate multiple components
	bonds = deepcopy(parent_substructure.bonds)
	for e in edges(bonds) # loop over bonds
		if get_prop(bonds, e, :cross_boundary) # find cross-PB bonds
			rem_edge!(bonds, e)
		end
	end
	# find conn comps
	conn_comps = connected_components(bonds)
	@info "Conn comps in bonds during conglomerate!" n=length(conn_comps)
	ref_comp_id = argmax(length.(conn_comps))
	# set reference atom
	ref_atom = conn_comps[ref_comp_id][1]
	# loop over non-reference components and rectify images
	for comp_id in 1:length(conn_comps)
		if comp_id == ref_comp_id
			continue
		end
		# get index of some atom in the non-reference comp
		atom_idx = conn_comps[comp_id][1]
		# find displacement vector for first atom
		dx = parent_substructure.atoms.coords.xf[:, ref_atom] - 
			 parent_substructure.atoms.coords.xf[:, atom_idx]
		# get nearest image vector
		nx = copy(dx)
		nearest_image!(nx)
		# if the norm of nx is less than that of dx, translate the atom
		@info "Vectors" dx nx
		if norm(nx) < norm(dx)
			@warn "Translating a component"
			for atom_idx in conn_comps[comp_id]
				parent_substructure.atoms.coords.xf[:, atom_idx] .+= dx - nx
			end
		end
	end
	return
end

# ╔═╡ 1ac0bd1b-5545-4ff8-bea4-aa88ca08281e
function get_r2p_alignment(replacement::Crystal, parent::Crystal, r2p::Dict{Int, Int})
    center = (X::Matrix{Float64}) -> sum(X, dims=2)[:] / size(X, 2)
    # when both centered to origin
    @assert replacement.atoms.n ≥ 3 && parent.atoms.n ≥ 3 "Parent and replacement must each be at least 3 atoms for SVD alignment."
    ###
    #   compute centered Cartesian coords of the atoms of 
    #       replacement fragment involved in alignment
    ###
    atoms_r = Cart(replacement.atoms[[r for (r, p) in r2p]], replacement.box)
    X_r = atoms_r.coords.x
    x_r_center = center(X_r)
    X_r = X_r .- x_r_center

    ###
    #   compute centered Cartesian coords of the atoms of 
    #       parent involved in alignment
    ###
    parent_substructure = deepcopy(parent[[p for (r, p) in r2p]])
    conglomerate!(parent_substructure)
    atoms_p = Cart(parent_substructure.atoms, parent_substructure.box)
    Xtals.write_xyz(atoms_p, "atoms_p.xyz")
    X_p = atoms_p.coords.x
    x_p_center = center(X_p)
    X_p = X_p .- x_p_center

    # solve the orthogonal procrustes probelm via SVD
    F = svd(X_r * X_p')
    # optimal rotation matrix
    rot =  F.V * F.U'

    err = norm(rot * X_r - X_p)

    return Alignment(rot, - x_r_center, x_p_center, err)
end

# ╔═╡ 66ee03a5-1635-4100-a5ad-53483cb4135e
md"## align replacement"

# ╔═╡ cc9d44ed-2b90-48df-83e9-7181ca2f8286
# do before calling optimal_replacement
begin
	nb_not_masked = sum(.! occursin.(rc[:r_tag], String.(query.atoms.species)))
	q_unmasked_in_r = query[1:nb_not_masked] ∈ replacement
	q2r = Dict([q => q_unmasked_in_r.isomorphisms[1][1][q] for q = 1:nb_not_masked])
end

# ╔═╡ 0dfe1c41-cadd-4d56-bdf2-106ac17e8df4
## these can probably be removed in the production version for efficiency
begin
	local is_masked = occursin.(rc[:r_tag], String.(query.atoms.species))
	@assert sum(is_masked) + nb_not_masked == query.atoms.n
	@assert ! any(is_masked[1:nb_not_masked]) "all masked atoms should be at the end"
end

# ╔═╡ 37a4feae-b957-44a5-8631-3811eecfaef5
md"## add crystals together"

# ╔═╡ a1cf7047-22dd-4747-82b9-2f1acfca730d
struct Installation
	aligned_replacement::Crystal
	q2p::Dict{Int, Int}
	r2p::Dict{Int, Int}
end

# ╔═╡ 1871494b-13cd-4845-9482-6f4afb4a62e7
function optimal_replacement(search::Search, replacement::Crystal, q2r::Dict{Int,Int}, loc_id::Int, ori_ids::Vector{Int}=[1:nb_ori_at_loc(search)[loc_id]...])
	# unpack search arg
	isomorphisms, parent, query = search.isomorphisms, search.parent, search.query

	# loop over ori_ids to find best r2p_alignment
	r2p_alignment = Alignment(zeros(1,1), [0.], [0.], Inf)
	best_ori = 0
	best_r2p = Dict{Int, Int}()
	for ori_id in ori_ids
		# find r2p isom
		q2p = isomorphisms[loc_id][ori_id]
		r2p = Dict([r => q2p[q] for (q, r) in q2r])
		# calculate alignment
		test_alignment = get_r2p_alignment(replacement, parent, r2p)
		# keep best alignment and generating ori_id
		if test_alignment.err < r2p_alignment.err
			r2p_alignment = test_alignment
			best_ori = ori_id
			best_r2p = r2p
		end
	end

	opt_aligned_replacement = aligned_replacement(replacement, parent, r2p_alignment)
	
	# return the replacement modified according to r2p_alignment
	return Installation(opt_aligned_replacement, isomorphisms[loc_id][best_ori], best_r2p)
end

# ╔═╡ fe0ec889-fca4-4d88-ba33-f800f73cc023
installation = optimal_replacement(search, replacement, q2r, 1)

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
	@info "keep_atoms" length(keep_atoms) child.atoms.n obsolete_atoms 26 in obsolete_atoms
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

# ╔═╡ 832d0c11-b3ec-498b-b2b1-cd2733c273de
parent.bonds

# ╔═╡ 0daa165f-dc4a-414f-bfb8-9a3716dca8c8
parent[1:end].bonds

# ╔═╡ 67aeb5ab-9dd6-481b-b1d5-ed214e3acfdc
+(parent, child, check_overlap=false).atoms.n

# ╔═╡ 647f4fcc-c5fb-4571-a5ab-daff9f10fef4
+(child, parent, check_overlap=false).atoms.n

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
# ╟─67692b25-bc22-4663-be60-8f3a5ee4a4cd
# ╠═1ac0bd1b-5545-4ff8-bea4-aa88ca08281e
# ╠═b289f607-6b82-4320-b326-533b982ec6c6
# ╟─66ee03a5-1635-4100-a5ad-53483cb4135e
# ╠═0dfe1c41-cadd-4d56-bdf2-106ac17e8df4
# ╠═cc9d44ed-2b90-48df-83e9-7181ca2f8286
# ╠═1871494b-13cd-4845-9482-6f4afb4a62e7
# ╠═fe0ec889-fca4-4d88-ba33-f800f73cc023
# ╟─37a4feae-b957-44a5-8631-3811eecfaef5
# ╠═a1cf7047-22dd-4747-82b9-2f1acfca730d
# ╠═94cfc1d5-5ee9-43d5-b936-0541d260f511
# ╠═586f38c8-4bd5-4cd5-a748-896f795afc6f
# ╟─f69bf2e7-87f8-403a-ae23-cf4b7bf7d9f1
# ╠═75e076ad-4cc0-42c3-9f2e-8138a810ae98
# ╠═832d0c11-b3ec-498b-b2b1-cd2733c273de
# ╠═0daa165f-dc4a-414f-bfb8-9a3716dca8c8
# ╠═67aeb5ab-9dd6-481b-b1d5-ed214e3acfdc
# ╠═647f4fcc-c5fb-4571-a5ab-daff9f10fef4
# ╠═9976fb79-210d-4eac-a99a-b81e02586039
# ╠═d8ea8d5e-f17a-412b-8461-15ba6d9621ec
# ╠═eae2225c-40f0-4d68-a9a2-43a39a82f029
