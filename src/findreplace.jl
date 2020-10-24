# findreplace.jl

## exposed interface
export substructure_search, SearchResult, Query, Search,
	nb_isomorphisms, nb_locations, nb_configs_at_loc, find_replace


## imports for extension
import Base.(∈)
import PorousMaterials.write_xyz


## Structs
# The parent structure and search moiety Crystals
struct Query
    parent::Crystal
    s_moty::Crystal
end

# The result of running a substructure search
struct Search
    query::Query # the search query (s-moty ∈ parent)
    results # the isomorphisms
end


## Helpers
@doc raw"""
`nb_isomorphisms(search::Search) -> Int`

Returns the number of isomorphisms found in the specified `Search`
"""
function nb_isomorphisms(search::Search)::Int
    return length(search.results.isomorphism)
end


@doc raw"""
`nb_locations(search::Search) -> Int`

Returns the number of unique locations (collections of atoms) at which the
specified `Search` results contain isomorphisms.
"""
function nb_locations(search::Search)::Int
    return length(unique(search.results.p_subset))
end


@doc raw"""
`nb_configs_at_loc(search::Search) -> Array{Int}`

Returns a array containing the number of isomorphic configurations at a given
location (collection of atoms) for which the specified `Search` results
contain isomorphisms.
"""
function nb_configs_at_loc(search::Search)::Array{Int}
    return [length(location.isomorphism) for location in groupby(search.results, :p_subset)]
end


# Retuns the geometric center of an Array, Frac/Atoms object, or Crystal.
function geometric_center(xf::Array{Float64,2})::Array{Float64}
	return sum(xf, dims=2)[:] / size(xf, 2)
end

geometric_center(coords::Frac)::Array{Float64} = geometric_center(coords.xf)

geometric_center(atoms::Atoms)::Array{Float64} = geometric_center(atoms.coords)

geometric_center(xtal::Crystal)::Array{Float64} = geometric_center(xtal.atoms)


# extension of infix `in` operator for expressive searching
# this allows all of the following:
#	s ∈ g					→	find the moiety in the crystal
#	[s1, s2] .∈ [g]			→	find each moiety in a crystal
#	s .∈ [g1, g2]			→	find the moiety in each crystal
#	[s1, s2] .∈ [g1, g2]	→	find each moiety in each crystal
(∈)(s::Crystal, g::Crystal) = substructure_search(s, g)
# and some more like that for doing find-replace operations in one line
# these are objectively unnecessary, but fun.
(∈)(pair::Pair{Crystal, Crystal}, xtal::Crystal) =
	find_replace((pair[1] ∈ xtal), pair[2], rand_all=true)
(∈)(tuple::Tuple{Pair, Int}, xtal::Crystal) =
	find_replace(tuple[1][1] ∈ xtal, tuple[1][2], nb_loc=tuple[2])
(∈)(tuple::Tuple{Pair, Array{Int}}, xtal::Crystal) =
	find_replace(tuple[1][1] ∈ xtal, tuple[1][2], loc=tuple[2])
(∈)(tuple::Tuple{Pair, Array{Int}, Array{Int}}, xtal::Crystal) =
	find_replace(tuple[1][1] ∈ xtal, tuple[1][2], loc=tuple[2], ori=tuple[3])


# Helper for making .xyz's
write_xyz(xtal::Crystal, name::String) = write_xyz(Cart(xtal.atoms, xtal.box), name)


# Translates all atoms in xtal such that xtal[1] is in its original position
# and the rest of xtal is in its nearest-image position relative to xtal[1]
function adjust_for_pb!(xtal::Crystal)
	# record position vector of xtal[1]
	origin_offset = deepcopy(xtal.atoms.coords.xf[:, 1])
	# loop over atom indices and correct coordinates
    for i in 1:xtal.atoms.n
		# move atoms near the origin for nearest-image calculation
		dxf = xtal.atoms.coords.xf[:, i] .- origin_offset
		# nearest_image! expects points to be within same or adjacent unit cells
		@assert all(abs.(dxf) .< 2) "Invalid xf coords in $(xtal.name)"
		# resolve periodic boundaries (other vectors unchanged)
        nearest_image!(dxf)
		# return atoms to their [nearest-image] original positions
        xtal.atoms.coords.xf[:, i] = dxf .+ origin_offset
    end
end


# Performs orthogonal Procrustes on correlated point clouds A and B
function orthogonal_procrustes(A::Array{Float64,2},
		B::Array{Float64,2})::Array{Float64,2}
	# solve the SVD
    F = svd(A * B')
	# return rotation matrix
    return F.V * F.U'
end


# Gets the s_moty-to-xtal rotation matrix
function s2p_op(s_moty::Crystal, xtal::Crystal)::Array{Float64,2}
	# s_moty in Cartesian
	A = s_moty.box.f_to_c * s_moty.atoms.coords.xf
	# parent subset in Cartesian
	B = xtal.box.f_to_c * xtal.atoms.coords.xf
	# get rotation matrix
	return orthogonal_procrustes(A, B)
end


# Gets the r_moty-to-s_mask rotation matrix
function r2m_op(r_moty::Crystal, s_moty::Crystal, m2r_isomorphism::Array{Int},
		s_mask_atoms::Atoms)::Array{Float64,2}
	# r_moty subset in Cartesian
	A = r_moty.box.f_to_c * r_moty.atoms[m2r_isomorphism].coords.xf
	# s_mask in Cartesian
	B = s_moty.box.f_to_c * s_mask_atoms.coords.xf
	# get rotation matrix
	return orthogonal_procrustes(A, B)
end


# Transforms r_moty according to two rotation matrices and a translational offset
function xform_r_moty(r_moty::Crystal, rot_r2m::Array{Float64,2},
		rot_s2p::Array{Float64,2}, xtal_offset::Array{Float64},
		xtal::Crystal)::Crystal
	# put r_moty into cartesian space
	atoms = Atoms{Cart}(length(r_moty.atoms.species), r_moty.atoms.species,
		Cart(r_moty.atoms.coords, r_moty.box))
	# transformation 1: rotate r_moty to align with s_moty
	atoms.coords.x[:,:] = rot_r2m * atoms.coords.x
	# transformation 2: rotate to align with xtal_subset
	atoms.coords.x[:,:] = rot_s2p * atoms.coords.x
	# transformation 3: translate to align with original xtal center
	atoms.coords.x .+= xtal.box.f_to_c * xtal_offset
	# cast atoms back to Frac
	return Crystal(r_moty.name, xtal.box, Atoms{Frac}(length(atoms.species),
		atoms.species, Frac(atoms.coords, xtal.box)), Charges{Frac}(0))
end


# shifts coordinates to make the geometric center of the point cloud coincident
# w/ the origin
function center_on_self!(xtal::Crystal)
	xtal.atoms.coords.xf .-= geometric_center(xtal)
end


# returns an Array containing the indices
function idx_filter(xtal::Crystal, subset::Array{Int})::Array{Int,1}
	return [i for i in 1:xtal.atoms.n if !(i ∈ subset)]
end


# tracks which bonds need to be made between the parent and array
# of transformed r_moty's (xrms)
function accumulate_bonds!(bonds::Array{Tuple{Int,Int}}, s2p_isom::Array{Int},
		parent::Crystal, m2r_isom::Array{Int}, r_moty::Crystal, xrms::Array{Crystal})
    # loop over s2p_isom
    for (s, p) in enumerate(s2p_isom)
        # find neighbors of parent_subset atoms
        n = LightGraphs.neighbors(parent.bonds, p)
        # loop over neighbors
        for nᵢ in n
            # if neighbor not in s2p_isom, must bond it to r_moty replacement
			# of parent_subset atom
            if !(nᵢ ∈ s2p_isom)
                # ID the atom in r_moty
                r = m2r_isom[s]

                # add the index offset
                r += parent.atoms.n + length(xrms) * r_moty.atoms.n

                # push bond to array
                push!(bonds, (nᵢ, r))
            end
        end
    end
end


# generates data for effecting a series of replacements
function build_replacement_data(c::Array{Int}, search::Search, parent::Crystal,
		s_moty::Crystal, r_moty::Crystal, m2r_isom::Array{Int}, mask::Crystal
        )::Tuple{Array{Crystal},Array{Int},Array{Tuple{Int,Int}}}
    xrms = Crystal[]
    del_atoms = Int[]
    bonds = Tuple{Int,Int}[] # each tuple (i,j) encodes a parent[i] -> xrms[k][j] bond
    for cᵢ in c
        # find isomorphism
        s2p_isom = search.results.isomorphism[cᵢ]
        # find parent subset
        parent_subset = deepcopy(parent[s2p_isom])
        # adjust coordinates for periodic boundaries
        adjust_for_pb!(parent_subset)
        # record the center of xtal_subset so we can translate back later
        parent_subset_center = geometric_center(parent_subset)
        # shift to align centers at origin
        center_on_self!.([parent_subset, s_moty])
        # do orthog. Procrustes for s_moty-to-parent and mask-to-replacement alignments
        rot_s2p = s2p_op(s_moty, parent_subset)
        rot_r2m = r2m_op(r_moty, s_moty, m2r_isom, mask.atoms)
        # transform r_moty by rot_r2m, rot_s2p, and xtal_subset_center, align to parent
        # (this is now a crystal to add)
        xrm = xform_r_moty(r_moty, rot_r2m, rot_s2p, parent_subset_center, parent)
        push!(xrms, xrm)
        # push obsolete atoms to array
        for x in s2p_isom
            push!(del_atoms, x) # this can create duplicates; remove them later
        end
        # clean up del_atoms
        del_atoms = unique(del_atoms)
        # accumulate bonds
        accumulate_bonds!(bonds, s2p_isom, parent, m2r_isom, r_moty, xrms)
    end
    return xrms, del_atoms, bonds
end


# returns ids of results corresponding to the
function get_config_ids(search::Search; orientations::Array{Int}=Int[], nb_loc::Int=0,
		locations::Array{Int}=Int[])::Array{Int}
    # randomize orientations
    if orientations == []
		orientations = [rand(1:n) for n ∈ nb_configs_at_loc(search)]
	end
    nb_loc = nb_loc == 0 ? nb_locations(search) : nb_loc
	if locations == []
		locations = [1:nb_locations(search)...]
	else
		nb_loc = length(locations)
	end
    # get id of ith isom in each group per orientations array
	ids = Int[]
	for (i, locᵢ) ∈ enumerate(groupby(search.results, :p_subset))
		if i ∈ locations && length(ids) < nb_loc
			push!(ids, locᵢ.id[orientations[length(ids)+1]])
		end
	end
	return ids
end


## Search function (exposed)
@doc raw"""
`substructure_search(s_moty::Crystal, xtal::Crystal) -> Search`

Searches for a substructure within a `Crystal` and returns a
`Search` struct containing all identified subgraph isomorphisms.
"""
function substructure_search(s_moty::Crystal, xtal::Crystal)::Search
	# Make a copy w/o R tags for searching
	moty = deepcopy(s_moty)
	untag_r_group!(moty)
	# Get array of configuration arrays
	configurations = find_subgraph_isomorphisms(moty.bonds,
		moty.atoms.species, xtal.bonds,	xtal.atoms.species)
    return Search(
		Query(xtal, s_moty),
		DataFrame(
			id = [1:length(configurations)...],
			p_subset = [sort(c) for c in configurations],
			isomorphism = configurations))
end


## Internal method for performing substructure replacements
function substructure_replace(s_moty::Crystal, r_moty::Crystal, parent::Crystal,
		search::Search, configs::Array{Int})::Crystal
	# configs must all be unique
	@assert length(configs) == length(unique(configs))
	# mutation guard
    s_moty, r_moty = deepcopy.([s_moty, r_moty])
    # if there are no replacements to be made, just return the parent
    if length(search.results.isomorphism) == 0
		@warn "No replacements to be made."
        return parent
    end
    # determine s_mask (which atoms from s_moty are NOT in r_moty?)
    mask = s_moty[idx_filter(s_moty, r_group_indices(s_moty))]
    # get isomrphism between s_moty/mask and r_moty
    m2r_isom = (mask ∈ r_moty).results.isomorphism[1]
    # shift all r_moty nodes according to center of isomorphic subset
    r_moty.atoms.coords.xf .-= geometric_center(r_moty[m2r_isom])
    # loop over configs to build replacement data
    xrms, del_atoms, bonds = build_replacement_data(configs, search, parent, s_moty,
		r_moty, m2r_isom, mask)

    # append temporary crystals to parent
    xtal = Crystal("TODO", parent.box,
		parent.atoms + sum([xrm.atoms for xrm in xrms]), Charges{Frac}(0))

    # create bonds from dictionary
    for (p, r) in bonds
        add_edge!(xtal.bonds, p, r)
    end
    # correct for periodic boundaries
	wrap!(xtal)
    # delete atoms from array and return result
    return xtal[[x for x in 1:xtal.atoms.n if !(x ∈ del_atoms)]]
end


## Find/replace function (exposed)
@doc raw"""
`find_replace(search::Search, r_moty::Crystal; rand_all::Bool=false, nb_loc::Int=0,
	loc::Array{Int}=Int[], ori::Array{Int}=Int[]) -> Crystal`

Inserts `r_moty` into a parent structure according to `search` and `kwargs`
"""
function find_replace(search::Search, r_moty::Crystal; rand_all::Bool=false,
		nb_loc::Int=0, loc::Array{Int}=Int[], ori::Array{Int}=Int[])::Crystal
    if rand_all
        return substructure_replace(search.query.s_moty, r_moty, search.query.parent,
			search, get_config_ids(search))
    elseif nb_loc > 0
        return substructure_replace(search.query.s_moty, r_moty, search.query.parent,
			search, get_config_ids(search, nb_loc=nb_loc));
	elseif ori ≠ Int[]
        return substructure_replace(search.query.s_moty, r_moty, search.query.parent,
			search, get_config_ids(search, locations=loc, orientations=ori))
	elseif loc ≠ Int[]
        return substructure_replace(search.query.s_moty, r_moty, search.query.parent,
			search, get_config_ids(search, locations=loc))
    else
        @error "No replacements specified."
    end
    return search.query.parent
end
