# findreplace.jl

## exposed interface

export substructure_search, find_replace, SearchResult, SearchQuery,
	Search, Configuration, convert, isequal


## imports for extension

import Base.convert
import Base.isequal
import Base.(∈)
import PorousMaterials.write_xyz


## Structs, conversions, and helpers

# for encoding the location and orientation indices within a group of search results
struct Configuration
    location::Int
    orientation::Int
end

# default case (invalid configuration)
Configuration() = Configuration(0, 0)

# hacks for weird bugs and convenient syntax
convert(Configuration, (i,j)) = Configuration(i,j)
Configuration(t::Tuple{Int,Int}) = Configuration(t[1],t[2])

# make Configurations display nicely BUG does not play nice with @info in jupyter
function Base.show(io::IO, c::Configuration)
	print("($(c.location),$(c.orientation))")
end

# for describing the search query ("find x in y")
struct SearchQuery
    substructure::Crystal
    parent_structure::Crystal
end

# empty default constructor
SearchQuery() = SearchQuery(nothing, nothing)

# hack for convenience
convert(SearchQuery,
	(a,b)::Tuple{Crystal, Crystal}) = SearchQuery(a,b)

# hack to display queries nicely
function Base.show(io::IO, sq::SearchQuery)
	print("$(sq.substructure.name) ∈ $(sq.parent_structure.name)")
end

# for storing search results in a user-friendly way
struct SearchResult
    query::SearchQuery
    configuration::Configuration
    isomorphism::Array{Int}
end

# make SearchResults look better
function Base.show(io::IO, ssr::SearchResult)
	print(ssr.query)
	print(" ")
	print(ssr.configuration)
	print(" ")
	print(ssr.isomorphism)
end

# stores a set of search results and basic meta data
struct Search
    results::Array{SearchResult}
    num_isomorphisms::Int
    num_locations::Int
    num_orientations::Array{Int}
end

# default constructor
Search() = Search([], 0, 0, [])

# display helper
function Base.show(io::IO, ssrs::Search)
	print("$(ssrs.num_isomorphisms) result(s) in $(ssrs.num_locations) location(s)")
	for result in ssrs.results
		print("\n\t")
		print(result)
	end
	print("\n")
end

# define equality for Search object
isequal(a::Search, b:: Search) =
	a.num_isomorphisms == b.num_isomorphisms &&
	a.num_locations == b.num_locations &&
	all(a.results .== b.results) &&
	all(a.num_orientations .== b.num_orientations)

## Helpers

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

function config_idx_arr(configs::Array{Configuration}, search::Search)::Array{Int}
    c = Int[]
    if configs == Configuration[] # if configs not specified, choose one randomly
        push!(c, rand(1:search.num_isomorphisms))
    else # find indices of specified configs
        for config in configs
            for i in 1:search.num_isomorphisms
                if search.results[i].configuration == config
                    push!(c, i)
                    break
                end
            end
        end
    end
    return c
end


function accumulate_bonds!(bonds::Array{Tuple{Int,Int}}, s2p_isom::Array{Int}, parent::Crystal, m2r_isom::Array{Int}, r_moty::Crystal, xrms::Array{Crystal})
    # loop over s2p_isom
    for (s, p) in enumerate(s2p_isom)
        # find neighbors of parent_subset atoms
        n = LightGraphs.neighbors(parent.bonds, p)
        # loop over neighbors
        for nᵢ in n
            # if neighbor not in s2p_isom, must bond it to r_moty replacement of parent_subset atom
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


function build_replacement_data(c::Array{Int}, search::Search, parent::Crystal, s_moty::Crystal, r_moty::Crystal, m2r_isom::Array{Int}, mask::Crystal
        )::Tuple{Array{Crystal},Array{Int},Array{Tuple{Int,Int}}}
    xrms = Crystal[]
    del_atoms = Int[]
    bonds = Tuple{Int,Int}[] # each tuple (i,j) encodes a parent[i] -> xrms[k][j] bond
    for cᵢ in c
        # find isomorphism
        s2p_isom = search.results[cᵢ].isomorphism
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


## Search function (exposed)

@doc raw"""
Searches for a substructure within a `Crystal` and returns a
`Search` struct containing all identified subgraph
isomorphisms.
"""
function substructure_search(s_moty::Crystal, xtal::Crystal)::Search
	# Make a copy w/o R tags for searching
	moty = deepcopy(s_moty)
	untag_r_group!(moty)
	# Get array of configuration arrays
	configurations = Ullmann.find_subgraph_isomorphisms(moty.bonds,
		moty.atoms.species, xtal.bonds,	xtal.atoms.species)
	# Group configurations by location: [node sets] .=> [configurations]
	results = Dict{Array, Array{Array{Int}}}()
	for configuration in configurations
		sorted_config = sort(configuration)
		if sorted_config ∈ keys(results)
			push!(results[sorted_config], configuration)
		else
			merge!(results, Dict([sorted_config => [configuration]]))
		end
	end
	# unpack dict
	num_orientations = [
		length(results[location_set]) for location_set in keys(results)]
	results = [result for result in values(results)]
	# repack as array of structs
	result_struct_array = Array{SearchResult}([])
    for (location, isomorphisms) ∈ enumerate(results)
        for (orientation, isomorphism) ∈ enumerate(results[location])
            push!(result_struct_array, SearchResult((s_moty, xtal),
				(location, orientation), isomorphism))
        end
    end
    return Search(result_struct_array,
		length(result_struct_array), length(keys(results)), num_orientations)
end


## Find/replace function (exposed)

@doc raw"""
Finds the search moiety in the parent structure and replaces it at a location.

			TODO: docs
"""
function find_replace(s_moty::Crystal, r_moty::Crystal, parent::Crystal;
        configs::Array{Configuration}=Configuration[],
        search::Search=Search(),
        kwargs...)::Crystal
    # mutation guard
    s_moty, r_moty = deepcopy.([s_moty, r_moty])
    # if search not specified, do it now
	search = isequal(search, Search()) ? s_moty ∈ parent : search
    # if there are no replacements to be made, just return the parent
    if search.num_isomorphisms == 0
        return parent
    end
    # determine configuration index array
    c = config_idx_arr(configs, search)
    # determine s_mask (which atoms from s_moty are NOT in r_moty?)
    mask = s_moty[idx_filter(s_moty, r_group_indices(s_moty))]
    # get isomrphism between s_moty/mask and r_moty
    m2r_isom = (mask ∈ r_moty).results[1].isomorphism
    # shift all r_moty nodes according to center of isomorphic subset
    r_moty.atoms.coords.xf .-= geometric_center(r_moty[m2r_isom])
    # loop over configs to build replacement data
    xrms, del_atoms, bonds = build_replacement_data(c, search, parent, s_moty, r_moty, m2r_isom, mask)
    # append temporary crystals to parent
    xtal = Crystal("TODO", parent.box, parent.atoms + sum([xrm.atoms for xrm in xrms]), Charges{Frac}(0))
    # create bonds from dictionary
    for (p, r) in bonds
        add_edge!(xtal.bonds, p, r)
    end
	# correct for periodic boundaries
	wrap!(xtal)
    # delete atoms from array and return result
    return xtal[[x for x in 1:xtal.atoms.n if !(x ∈ del_atoms)]]
end
