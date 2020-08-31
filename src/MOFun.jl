module MOFun

## packages that are used by this module
using PorousMaterials
using LightGraphs
using LinearAlgebra
using Printf
using DataFrames
using Ullmann
include("moiety.jl")
## other necessary files containing function definitions and tests
include("ring_constructor.jl")
include("alignment_operations.jl")

## file paths for fragment files
fragment_location = joinpath(pwd(), "fragments")

## function that will make it easier to name files
remove_extension(crystal::Crystal) = split(crystal.name, ".")[1]

"""
Create a Crystal from a .xyz file with box=unit_cube()
and charges = Charges(0, Float64[], Frac(Array{Float64}(undef,3,0))).
"""
function read_fragment_from_xyz(name; fragment_location=fragment_location)
    filename = joinpath(fragment_location, name * ".xyz")
    atoms_c = read_xyz(filename)
    # contains: number of atoms, atom type (symbol), and position in fraction coords
    atoms_f = Frac(atoms_c, unit_cube())
    charges_f = Charges{Frac}(0) # create an empty charge struct
    return Crystal(name, unit_cube(), atoms_f, charges_f)
end

# The side of the ring wish to functionalize
# must be an Int (if substitution_position == "all" this parameter won't matter)
function choose_side(randomize_side::Bool=true, set_side::Int=2)
    # The side of the ring wish to functionalize
    # must be an Int (if substitution_position == "all" this parameter won't matter)
    if randomize_side
        return which_side = rand((1, 2))
    elseif !randomize_side
        return which_side = set_side
    end
end


## main function
function functionalize_mof(crystal::Crystal, fragment_name::String, ipso_species::Symbol, r_species::Symbol,
				bonding_rules::Array{BondingRule,1}; n::Int=6,
				side_to_functionalize::Int=2, randomize_side::Bool=true,
				arene_substitution_type::String="meta")
	## check for/create a directory to store the output files
	if ! isdir(joinpath(remove_extension(crystal)))
		mkdir(joinpath(remove_extension(crystal)))
	end

	####
	# initialize parameters
	####
	# difine the fragment to be used for the functionalization
	fragment = read_fragment_from_xyz(fragment_name)
	infer_bonds!(fragment, false, bonding_rules)

	# define what type of arene substitution to perform
	substitution_position = arene_substitution_type

	####
	# create an array populated with AromaticRings from the Crystal
	####
	rings = find_rings(crystal, r_species, ipso_species)

	####
	# perform alignment
	####
	aligned_fragments = []
	fragment_r_ids = []
	atoms_to_remove = Array{Int, 1}(undef, 0)
	fragment_atoms_to_remove = Array{Int, 1}(undef, 0)
	bonds_to_make = Array{Array{Int,1}, 1}(undef, 0)

	for (ring_id, ring) in enumerate(rings)
		# determine which side of the ring to functionalize
		which_side = choose_side(randomize_side, side_to_functionalize)

		# align fragment with appropriate place on ring in MOF
    		aligned_frag = align_fragment(ring, crystal, fragment,
				substitution_position, which_side=which_side)
    		push!(aligned_fragments, aligned_frag)

    		# identify the corsponding fragment_r_id with adjusted index
    		pos = (ring_id - 1) * fragment.atoms.n # number to adjust the indexes by
    		push!(fragment_r_ids,
        		crystal.atoms.n + pos .+ fragment_aro_R_id(aligned_frag))
    		append!(fragment_atoms_to_remove,
        		crystal.atoms.n + pos .+ fragment_aro_triplet_ids(aligned_frag))

    		# identify the hydrogen from the crystal attached to the C_aro_r
		# (the C we will functionalize)
    		if substitution_position == "meta"
        		C_to_bond = ring.sides[which_side].meta.id_C
        		push!(atoms_to_remove, ring.sides[which_side].meta.id_H)
    		elseif substitution_position == "ortho"
        		C_to_bond = ring.sides[which_side].ortho.id_C
        		push!(atoms_to_remove, ring.sides[which_side].ortho.id_H)
    		end
    		push!(bonds_to_make, [C_to_bond, fragment_r_ids[ring_id]])
	end

	####
	# connect Fragment to the Crystal by adding appropriate bonds
	# then, remove appropriate atoms and their bonds
	####
	append!(atoms_to_remove, fragment_atoms_to_remove)

	all_frags = +(aligned_fragments..., check_overlap=false)

	dirty_functionalized_mof = +(crystal, all_frags, check_overlap=false)

	for ed in bonds_to_make
		add_edge!(dirty_functionalized_mof.bonds, ed[1], ed[2])
	end

	functionalized_mof = dirty_functionalized_mof[[i for i = 1:dirty_functionalized_mof.atoms.n
				if ! (i in atoms_to_remove)]]
	####
	# write output files
	####
	write_xyz(Cart(functionalized_mof.atoms, functionalized_mof.box),
		joinpath(remove_extension(crystal),
		remove_extension(crystal) * "_" * substitution_position * "_functionalized_" * fragment.name))

	write_bond_information(functionalized_mof,
		joinpath(remove_extension(crystal),
		remove_extension(crystal) * "_" * substitution_position * "_functionalized_" *
		fragment.name * "_bonds.vtk"))

	final_cif_name = joinpath(remove_extension(crystal),
	split(crystal.name, ".")[1] * "_" * substitution_position * "_functionalized_" *
	fragment.name * ".cif")

	functionalized_mof = apply_symmetry_operations(functionalized_mof)
	functionalized_mof.atoms.coords.xf .= mod.(functionalized_mof.atoms.coords.xf, 1.0)

	write_cif(functionalized_mof, final_cif_name)
end


# for encoding the location and orientation indices within a group of search results
struct Configuration
    location::Int
    orientation::Int
end
import Base.convert
convert(Configuration, (i,j)) = Configuration(i,j)


# for describing the search query ("find x in y")
struct SubstructureQuery
    substructure::String
    parent_structure::String
    # Should these be Crystals?
    # Pro: easier to do analysis later
    # Con: memory bloat
end
convert(SubstructureQuery, (a,b)::Tuple{Crystal, Crystal}) = SubstructureQuery(a.name,b.name)


# for storing search results in a user-friendly way
struct SubstructureSearchResult
    query::SubstructureQuery
    configuration::Configuration
    isomorphism::Array{Int}
end


# stores a set of search results and basic meta data
struct SubstructureSearchResults
    results::Array{SubstructureSearchResult}
    num_isomorphisms::Int
    num_locations::Int
    num_orientations::Array{Int}
end


@doc raw"""
Wraps find_subgraph_isomorphisms for convenient substructure searching
"""
function substructure_search(find_moiety::Crystal,
	parent_structure::Crystal)::SubstructureSearchResults
	# Make a copy w/o R tags for searching
	moty = deepcopy(find_moiety)
	untag_r_group!(moty)
	# Get array of configuration arrays
	configurations = Ullmann.find_subgraph_isomorphisms(moty.bonds,
											  moty.atoms.species,
											  parent_structure.bonds,
											  parent_structure.atoms.species)
	# Group configurations by location
	# keys are sets of nodes describing a location
	# values are collections of configuration arrays
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
	location_sets = keys(results)
	num_orientations = [length(results[location_set]) for location_set in location_sets]
	results = [result for result in values(results)]
	# repack as array of structs. refactor as list comp?
	result_struct_array = Array{SubstructureSearchResult}([])
    for (location, isomorphisms) ∈ enumerate(results)
        @debug "Location $location"
        for (orientation, isomorphism) ∈ enumerate(results[location])
            push!(result_struct_array, SubstructureSearchResult(
				(find_moiety, parent_structure),
				(location,orientation),isomorphism))
        end
    end
    return SubstructureSearchResults(result_struct_array,
		length(result_struct_array),length(location_sets),num_orientations)
end


# extension of infix `in` operator for expressive searching
import Base.∈
(∈)(s::Crystal, g::Crystal) = substructure_search(s, g)
# this allows all of the following:
#	s ∈ g					→	find the moiety in the crystal
#	[s1, s2] .∈ [g]			→	find each moiety in a crystal
#	s .∈ [g1, g2]			→	find the moiety in each crystal
#	[s1, s2] .∈ [g1, g2]	→	find each moiety in each crystal


function adjust_for_pb!(atoms_array::Array{Atoms{Frac}})
   for atoms in atoms_array
        atoms2 = deepcopy(atoms)
        for i in 1:length(atoms.species)
            dxf = atoms.coords.xf[:, i] .- atoms.coords.xf[:, 1] # rebuild the image (fixing first node @ origin)
            nearest_image!(dxf)
            atoms2.coords.xf[:, i] = dxf
        end
        atoms= atoms2
    end
end



geometric_center(coords) = sum(coords, dims=2)[:] / 3


function orthogonal_procrustes(A, B)
    F = svd(A * B')
    return F.V * F.U'
end



function perform_ops(s_moty, r_moty, xtal, s2p_isomorphism, r_moty_subset_coords, s_mask_coords)
    rot_s2p = orthogonal_procrustes(s_moty.box.f_to_c*s_moty.atoms.coords.xf,
        xtal.box.f_to_c*xtal.atoms.coords.xf[:,s2p_isomorphism])
    rot_r2m = orthogonal_procrustes(r_moty.box.f_to_c*r_moty_subset_coords,
        s_moty.box.f_to_c*s_mask_coords)
    return rot_s2p, rot_r2m
end



function shift_to_origin!(coord_arrays::Array{Array{Float64,2}})
    for array in coord_arrays
        array .-= geometric_center(array)
    end
end



function xform_r_moty(r_moty, rot_r2m, rot_s2p, xtal_subset_center, xtal)
    xformd_r_moty_coords = rot_r2m * r_moty.box.f_to_c * r_moty.atoms.coords.xf
    xformd_r_moty_atoms = Atoms{Cart}(length(r_moty.atoms.species),
        r_moty.atoms.species,
        Cart(rot_s2p * xformd_r_moty_coords .+ xtal.box.f_to_c * xtal_subset_center))
    return Frac(xformd_r_moty_atoms, xtal.box)
end



const DEBUG = true
function effect_replacement(xtal::Crystal, s_moty::Crystal, r_moty::Crystal, s2p_isomorphism::Array{Int}, m2r_isomorphism)::Crystal
    if DEBUG
        write_xyz(Cart(xtal.atoms, xtal.box), "01_xtal")
        write_xyz(Cart(r_moty.atoms, r_moty.box), "02_r_moty")
        write_xyz(Cart(s_moty.atoms, s_moty.box), "03_s_moty")
    end
    # adjust atomic coordinates to account for periodic boundaries
    xtal_subset = xtal.atoms[s2p_isomorphism] # adjust coords on this, to preserve rest of original crystal in unit cell
    xtal_subset_center = geometric_center(xtal_subset.coords.xf) # center of location for replacement
    adjust_for_pb!([xtal_subset, s_moty.atoms, r_moty.atoms])
    if DEBUG
        write_xyz(Cart(xtal_subset, xtal.box), "04_adjusted_xtal_subset")
        write_xyz(Cart(s_moty.atoms, s_moty.box), "05_adjusted_s_moty")
        write_xyz(Cart(r_moty.atoms, r_moty.box), "06_adjusted_r_moty")
    end

    # determine s_mask
    r_indices = MOFun.r_group_indices(s_moty) # which atoms from s_moty are NOT in r_moty?
    s_mask_indices = [index for index in 1:length(s_moty.atoms.species) if !(index ∈ r_indices)]
    s_mask_atoms = s_moty.atoms[s_mask_indices]
    s_mask_coords = s_mask_atoms.coords.xf
    if DEBUG
        write_xyz(Cart(s_moty.atoms[r_indices], s_moty.box), "07_r_group")
        write_xyz(Cart(s_mask_atoms, s_moty.box), "08_search_moiety_mask")
    end

    # shift coords to align centers at origin
    shift_to_origin!([s_moty.atoms.coords.xf, xtal_subset.coords.xf, s_mask_coords])
    # r_moty is different: shift all nodes according to center of a subset
    r_moty.atoms.coords.xf .-= geometric_center(r_moty.atoms.coords.xf[:, m2r_isomorphism])
    if DEBUG
        write_xyz(Cart(s_moty.atoms, s_moty.box), "09_shifted_s_moty")
        write_xyz(Cart(xtal_subset, xtal.box), "10_shifted_xtal_subset")
        write_xyz(Cart(r_moty.atoms, r_moty.box), "11_shifted_r_moty")
    end

    # do orthogonal Procrustes for s_moty-to-parent and mask-to-replacement alignments
    rot_s2p, rot_r2m = perform_ops(s_moty, r_moty, xtal, s2p_isomorphism,
        r_moty.atoms[m2r_isomorphism].coords.xf, s_mask_coords)

    # transform r_moty according to rot_r2m, rot_s2p, and xtal_subset_center, align to box
    r_coords = xform_r_moty(r_moty, rot_r2m, rot_s2p, xtal_subset_center, xtal)

    # subtract s_moty isomorphic subset from xtal
    keep = [i for i in 1:length(xtal.atoms.species) if !(i ∈ s2p_isomorphism)]
    crystal_species = xtal.atoms.species[keep]
    crystal_coords = xtal.atoms.coords[keep]

    # concatenate transformed r_moty's coords/species w/ xtal's
    species = vcat(crystal_species, r_moty.atoms.species)
    coords = hcat(crystal_coords.xf, r_coords.coords.xf)

    # return new crystal
    name = remove_extension(xtal) * "_find_" * MOFun.remove_path_prefix(s_moty.name) * "_replace_" * r_moty.name
    atoms = Atoms(species, Frac(coords))
    crystal = Crystal(name, xtal.box, atoms, Charges{Frac}(0))
    wrap!(crystal)
    infer_bonds!(crystal, true)
    return crystal
end


"""
Finds the search moiety in the parent structure and replaces it at one or more locations.
"""
function find_replace(s_moty::Crystal, r_moty::Crystal, parent::Crystal; c::Int=1)::Crystal
	s2p_isomorphism = (s_moty ∈ parent).results[c].isomorphism
	mask = subtract_r_group(s_moty)
	infer_bonds!(mask, false)
	m2r_isomorphism = (mask ∈ r_moty).results[1].isomorphism # may eventually want others, if they occur?
	return effect_replacement(parent, s_moty, r_moty, s2p_isomorphism, m2r_isomorphism)
end


export
	# MOFfun.jl
	remove_extension, read_fragment_from_xyz, functionalize_mof,
	choose_side, substructure_search, SubstructureSearchResult,
	SubstructureQuery, SubstructureSearchResults, Configuration,
	find_replace,

	# ring_constructor.jl
	empty_ring, is_aromatic, find_aromatic_cycles,
	is_species, have_neighbor, find_ipso_id,
	cycle_to_ring, find_rings,

	# alignment_operations.jl
	crystal_aro_triplet_ids, fragment_aro_triplet_ids, fragment_aro_R_id,
	center_of_aro_triplet, triplet_locality, align_fragment,

	# moiety.jl
	moiety
end
