### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ cdf8f316-54a1-11ec-0e6d-b1a485ca5fe8
import Pkg; Pkg.activate()

# ╔═╡ 1860c789-ca2a-4d2e-be9f-af91fe544691
using PoreMatMod, Bio3DView, PlutoUI, Graphs, LinearAlgebra, StatsBase

# ╔═╡ b6dea881-26f4-4a20-b916-026dcae240c3
function viz(xtal::Crystal)
	write_cif(xtal, "dude.cif")
	run(`avogadro dude.cif`)
end

# ╔═╡ 536e88c0-3924-4cfb-a4f9-e5ad8c76ea88
parent = replicate(Crystal("NiPyC_fragment_trouble.cif"), (2,2,2))

# ╔═╡ f313cbbc-b39f-4ba1-ba25-6b106a9c5c0e
infer_bonds!(parent, true)

# ╔═╡ 1362ef07-97e2-4d6d-8c6e-5fbb56e55dbc
query = moiety("PyC.xyz")

# ╔═╡ acf39a52-4a36-49e2-802a-5c3960180246
replacement = moiety("PyC-CH3.xyz")

# ╔═╡ cf65e65a-6d7a-45f1-bdb3-a08c77f00c80
search = query in parent

# ╔═╡ 4f8868ec-d566-4525-af7e-c54801e3d7b9
begin
	# arguments
	search.isomorphisms
	loc_id = 1
	ori_id = 1
end

# ╔═╡ 9585bd79-1fd6-4086-b4da-7e65d6400fea
begin
	# TODO make this a dict by default
	q2p = Dict([q => search.isomorphisms[loc_id][ori_id][q] for q = 1:length(search.isomorphisms[loc_id][ori_id])])
end

# ╔═╡ b65a6988-40e1-4311-ac82-c9d6f03096ae
function ids_not_masked(query::Crystal)
	return filter(i -> ! occursin(rc[:r_tag], String(query.atoms.species[i])), 1:query.atoms.n)
end

# ╔═╡ fe31629c-ce57-46ff-9f76-9b09b1257164
nb_masked(query::Crystal) = sum(occursin.(rc[:r_tag], String.(query.atoms.species)))

# ╔═╡ cc9d44ed-2b90-48df-83e9-7181ca2f8286
nb_not_masked(query::Crystal) = sum(.! occursin.(rc[:r_tag], String.(query.atoms.species)))

# ╔═╡ 0dfe1c41-cadd-4d56-bdf2-106ac17e8df4
@assert nb_masked(query) + nb_not_masked(query) == query.atoms.n

# ╔═╡ 3b81670a-b4d3-4d03-a42c-1d033edffb93
@assert ids_not_masked(query) == 1:query.atoms.n-nb_masked(query)

# ╔═╡ 21a25153-9ce8-419f-8c3d-69cdeb0a0d67
q_unmasked_in_r = query[ids_not_masked(query)] ∈ replacement

# ╔═╡ a47dc0bd-73a9-49da-9253-9fc28b78be9f
q2r = Dict([q => q_unmasked_in_r.isomorphisms[1][1][q] for q = 1:nb_not_masked(query)])

# ╔═╡ 96f583de-ff3b-4287-bf16-132eaf9adc05
r2p = Dict([r => q2p[q] for (q, r) in q2r])

# ╔═╡ b289f607-6b82-4320-b326-533b982ec6c6
function conglomerate!(parent_substructure::Crystal)
	ref_atom = 1
	# todo
	return
end

# ╔═╡ 92bc3f84-3fab-4fa9-aee9-e94cafe1d9be
center(X::Matrix{Float64}) = sum(X, dims=2)[:] / size(X, 2)

# ╔═╡ e561a0f6-eaad-4ac8-b814-c85ed3b92a09
struct Alignment
	rot::Matrix{Float64}
	# before rotation
	shift_1::Vector{Float64}
	# after rotation
	shift_2::Vector{Float64}
	# error
	err::Float64
end

# ╔═╡ ab6d2e0c-bf7b-4513-a476-f9e52e244d39
function get_r2p_alignment(replacement::Crystal, 
	                       parent::Crystal, 
	                       r2p::Dict{Int, Int})
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
	parent_substructure = parent[[p for (r, p) in r2p]]
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

# ╔═╡ 04e61f68-7336-4c1d-888b-50fba2cae00d
r2p_alignment = get_r2p_alignment(replacement, parent, r2p)

# ╔═╡ 22ae2296-fe71-4e55-9e02-8b951bc07f59
function aligned_replacement(replacement::Crystal, 
				             parent::Crystal,
	                         r2p_alignment::Alignment
)
    # put replacement into cartesian space
    atoms_r = Cart(replacement.atoms, replacement.box)
    # rotate replacement to align with parent_subset
    atoms_r.coords.x[:, :] = r2p_alignment.rot * (atoms_r.coords.x .+ r2p_alignment.shift_1) .+ r2p_alignment.shift_2
    # cast atoms back to Frac
    return Crystal(replacement.name, parent.box, 
		           Frac(atoms_r, parent.box), Charges{Frac}(0))
end

# ╔═╡ 0cbf3baa-24fc-4022-a08e-8b24379d1500
replacement_to_install = aligned_replacement(replacement, parent, r2p_alignment)

# ╔═╡ 2b342469-8967-4817-9439-ea6625321c41
trim(parent::Crystal, q2p::Dict{Int, Int}) = parent[[p for p = 1:parent.atoms.n if ! (p in values(q2p))]]

# ╔═╡ 586f38c8-4bd5-4cd5-a748-896f795afc6f
child = trim(parent, q2p) + replacement_to_install

# ╔═╡ d8ea8d5e-f17a-412b-8461-15ba6d9621ec
write_cif(child)

# ╔═╡ 6168b10b-40db-439e-86df-179ea2a1c506
viz(child)

# ╔═╡ Cell order:
# ╠═cdf8f316-54a1-11ec-0e6d-b1a485ca5fe8
# ╠═1860c789-ca2a-4d2e-be9f-af91fe544691
# ╠═b6dea881-26f4-4a20-b916-026dcae240c3
# ╠═536e88c0-3924-4cfb-a4f9-e5ad8c76ea88
# ╠═f313cbbc-b39f-4ba1-ba25-6b106a9c5c0e
# ╠═1362ef07-97e2-4d6d-8c6e-5fbb56e55dbc
# ╠═acf39a52-4a36-49e2-802a-5c3960180246
# ╠═cf65e65a-6d7a-45f1-bdb3-a08c77f00c80
# ╠═4f8868ec-d566-4525-af7e-c54801e3d7b9
# ╠═9585bd79-1fd6-4086-b4da-7e65d6400fea
# ╠═b65a6988-40e1-4311-ac82-c9d6f03096ae
# ╠═fe31629c-ce57-46ff-9f76-9b09b1257164
# ╠═cc9d44ed-2b90-48df-83e9-7181ca2f8286
# ╠═0dfe1c41-cadd-4d56-bdf2-106ac17e8df4
# ╠═3b81670a-b4d3-4d03-a42c-1d033edffb93
# ╠═21a25153-9ce8-419f-8c3d-69cdeb0a0d67
# ╠═a47dc0bd-73a9-49da-9253-9fc28b78be9f
# ╠═96f583de-ff3b-4287-bf16-132eaf9adc05
# ╠═b289f607-6b82-4320-b326-533b982ec6c6
# ╠═92bc3f84-3fab-4fa9-aee9-e94cafe1d9be
# ╠═e561a0f6-eaad-4ac8-b814-c85ed3b92a09
# ╠═ab6d2e0c-bf7b-4513-a476-f9e52e244d39
# ╠═04e61f68-7336-4c1d-888b-50fba2cae00d
# ╠═22ae2296-fe71-4e55-9e02-8b951bc07f59
# ╠═0cbf3baa-24fc-4022-a08e-8b24379d1500
# ╠═2b342469-8967-4817-9439-ea6625321c41
# ╠═586f38c8-4bd5-4cd5-a748-896f795afc6f
# ╠═d8ea8d5e-f17a-412b-8461-15ba6d9621ec
# ╠═6168b10b-40db-439e-86df-179ea2a1c506
