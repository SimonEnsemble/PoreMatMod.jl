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
convert(SubstructureQuery,
	(a,b)::Tuple{Crystal, Crystal}) = SubstructureQuery(a.name,b.name)


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
Searches for a substructure within a `Crystal` and returns a `SubstructureSearchResults`
struct containing all identified subgraph isomorphisms.
"""
function substructure_search(find_moiety::Crystal,
	parent_structure::Crystal)::SubstructureSearchResults
	# Make a copy w/o R tags for searching
	moty = deepcopy(find_moiety)
	untag_r_group!(moty)
	# Get array of configuration arrays
	configurations = Ullmann.find_subgraph_isomorphisms(moty.bonds,
		moty.atoms.species, parent_structure.bonds,	parent_structure.atoms.species)
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
				(location, orientation), isomorphism))
        end
    end
    return SubstructureSearchResults(result_struct_array,
		length(result_struct_array), length(location_sets), num_orientations)
end


# extension of infix `in` operator for expressive searching
import Base.∈
(∈)(s::Crystal, g::Crystal) = substructure_search(s, g)
# this allows all of the following:
#	s ∈ g					→	find the moiety in the crystal
#	[s1, s2] .∈ [g]			→	find each moiety in a crystal
#	s .∈ [g1, g2]			→	find the moiety in each crystal
#	[s1, s2] .∈ [g1, g2]	→	find each moiety in each crystal


# Translates all atoms in xtal such that Atoms[1] is @ origin, w/ correction for
# the Crystal's periodic boundaries.
function adjust_for_pb(xtal::Crystal)::Crystal
	@debug "adjust_for_pb" xtal.name
    atoms = deepcopy(xtal.atoms)
    for i in 1:atoms.n
		 # rebuild the image (fixing first node @ origin)
		dxf = atoms.coords.xf[:, i] .- xtal.atoms.coords.xf[:, 1]
		# nearest_image! expects points to be within same or adjacent unit cells
		@assert all(abs.(dxf) .< 2) "Invalid xf coords in $(xtal.name)"
        nearest_image!(dxf)
        atoms.coords.xf[:, i] = dxf
    end
    return Crystal(xtal.name, xtal.box, atoms, xtal.charges)
end


# Retuns the geometric center of an array of points.
geometric_center(xf::Array{Float64,2})::Array{Float64} = sum(xf, dims=2)[:] / 3


# Performs orthogonal Procrustes on correlated point clouds A and B
function orthogonal_procrustes(A::Array{Float64,2}, B::Array{Float64,2})::Array{Float64,2}
    F = svd(A * B')
    return F.V * F.U'
end


# Gets the s_moty-to-xtal rotation matrix
function s2p_op(s_moty::Crystal,
		xtal::Crystal,
		s2p_isomorphism::Array{Int})::Array{Float64}
	A = s_moty.box.f_to_c * s_moty.atoms.coords.xf
	B = xtal.box.f_to_c * xtal.atoms.coords.xf[:, s2p_isomorphism]
	return orthogonal_procrustes(A, B)
end


# Gets the r_moty-to-s_mask rotation matrix
function r2m_op(r_moty::Crystal,
		s_moty::Crystal,
		m2r_isomorphism::Array{Int},
		s_mask_atoms::Atoms)::Array{Float64}
	r_moty_subset_xf = r_moty.atoms[m2r_isomorphism].coords.xf
	A = r_moty.box.f_to_c * r_moty_subset_xf
	B = s_moty.box.f_to_c * s_mask_atoms.coords.xf
	return orthogonal_procrustes(A, B)
end


# Shifts an array of points such that they are centered about the origin
function shift_to_origin!(xf::Array{Float64,2})
    xf .-= geometric_center(xf)
end


# Transforms r_moty according to two rotation matrices and a translational offset
function xform_r_moty(r_moty::Crystal,
		rot_r2m::Array{Float64,2},
		rot_s2p::Array{Float64,2},
		xtal_subset_center::Array{Float64},
		xtal::Crystal)
	# put r_moty into cartesian space
	coords = Cart(r_moty.atoms.coords, r_moty.box)
	atoms = Atoms{Cart}(length(r_moty.atoms.species), r_moty.atoms.species, coords)
	# transformation 1: rotate r_moty to align with s_moty
	atoms.coords.x[:,:] = rot_r2m * atoms.coords.x
	if DEBUG
		write_xyz(atoms, "output/12_rotation1.xyz")
	end
	# transformation 2: rotate to align with xtal_subset
	atoms.coords.x[:,:] = rot_s2p * atoms.coords.x
	if DEBUG
		write_xyz(atoms, "output/13_rotation2.xyz")
	end
	# transformation 3: translate to align with xtal
	atoms.coords.x .+= xtal.box.f_to_c * xtal_subset_center
	if DEBUG
		write_xyz(atoms, "output/14_translation.xyz")
	end
	# cast atoms back to Frac
	atoms = Atoms{Frac}(length(atoms.species), atoms.species, Frac(atoms.coords, xtal.box))
	if DEBUG
		write_xyz(Cart(atoms, xtal.box), "output/15_fractional.xyz")
	end
	return atoms
end


# Helper for making debugging .xyz's
import PorousMaterials.write_xyz
write_xyz(xtal::Crystal, name::String) = write_xyz(Cart(xtal.atoms, xtal.box), name)


# returns a Crystal with a specific replacement made
const DEBUG = true
function effect_replacement(xtal::Crystal,
		s_moty::Crystal,
		r_moty::Crystal,
		s2p_isomorphism::Array{Int},
		m2r_isomorphism)::Crystal
	# prevent input mutation
	s_moty = deepcopy(s_moty)
	r_moty = deepcopy(r_moty)
    if DEBUG # export inputs
        write_xyz(xtal, "output/01_xtal")
        write_xyz(r_moty, "output/02_r_moty")
        write_xyz(s_moty, "output/03_s_moty")
    end
    # adjust atomic coordinates to account for periodic boundaries
	# adjust coords on this, to preserve rest of original crystal in unit cell
    xtal_subset = deepcopy(xtal[s2p_isomorphism])
	# center of location for replacement
    xtal_subset_center = geometric_center(xtal_subset.atoms.coords.xf)
    xtal_subset, s_moty, r_moty = adjust_for_pb.([xtal_subset, s_moty, r_moty])
    if DEBUG # export inputs adjusted for periodic boundaries (atom 1 @ origin)
        write_xyz(xtal_subset, "output/04_adjusted_xtal_subset")
        write_xyz(s_moty, "output/05_adjusted_s_moty")
        write_xyz(r_moty, "output/06_adjusted_r_moty")
    end
## BUG moiety inputs get mangled by nearest_image! b/c they are not in Frac([0,1]³)

    # determine s_mask
    r_indices = r_group_indices(s_moty) # which atoms from s_moty are NOT in r_moty?
    s_mask_indices = [index for index in 1:length(s_moty.atoms.species)
		if !(index ∈ r_indices)]
    s_mask_atoms = s_moty.atoms[s_mask_indices]
    if DEBUG
        write_xyz(s_moty[r_indices], "output/07_r_group")
        write_xyz(Cart(s_mask_atoms, s_moty.box), "output/08_search_moiety_mask")
    end

    # shift coords to align centers at origin
    shift_to_origin!.([s_moty.atoms.coords.xf, xtal_subset.atoms.coords.xf,
		s_mask_atoms.coords.xf])
    # r_moty is different: shift all nodes according to center of a subset
    r_moty.atoms.coords.xf .-= geometric_center(r_moty.atoms.coords.xf[:, m2r_isomorphism])
    if DEBUG
        write_xyz(s_moty, "output/09_shifted_s_moty")
        write_xyz(xtal_subset, "output/10_shifted_xtal_subset")
        write_xyz(r_moty, "output/11_shifted_r_moty")
    end

    # do orthogonal Procrustes for s_moty-to-parent and mask-to-replacement alignments
    rot_s2p = s2p_op(s_moty, xtal, s2p_isomorphism)
	rot_r2m = r2m_op(r_moty, s_moty, m2r_isomorphism, s_mask_atoms)

    # transform r_moty according to rot_r2m, rot_s2p, and xtal_subset_center, align to box
    atoms = xform_r_moty(r_moty, rot_r2m, rot_s2p, xtal_subset_center, xtal)
	r_moty = Crystal(r_moty.name, xtal.box, atoms, Charges{Frac}(0))
	if DEBUG
		write_xyz(r_moty, "output/17_transformed.xyz")
	end

    # subtract s_moty isomorphic subset from xtal
    keep = [i for i in 1:length(xtal.atoms.species) if !(i ∈ s2p_isomorphism)]

	# put xtal and r_moty atoms into Cart, add together
	xtal_atoms_x = Cart(xtal.atoms[keep], xtal.box)
	r_moty_atoms_x = Cart(r_moty.atoms, r_moty.box)
	atoms = xtal_atoms_x + r_moty_atoms_x
	# put atoms back into Frac
	atoms = Frac(atoms, xtal.box)
	name = remove_extension(xtal) * "_find_" * MOFun.remove_path_prefix(s_moty.name) *
		"_replace_" * r_moty.name
    crystal = Crystal(name, xtal.box, atoms, Charges{Frac}(0))
    wrap!(crystal)
    infer_bonds!(crystal, true) # this needs to be changed to intelligent bond creation!
    return crystal
end


@doc raw"""
Finds the search moiety in the parent structure and replaces it at one or more locations.
"""
function find_replace(s_moty::Crystal, r_moty::Crystal, parent::Crystal; c::Int=1)::Crystal
	s2p_isomorphism = (s_moty ∈ parent).results[c].isomorphism
	mask = subtract_r_group(s_moty)
	infer_bonds!(mask, false)
	# may eventually want others, if they occur?
	m2r_isomorphism = (mask ∈ r_moty).results[1].isomorphism
	return effect_replacement(parent, s_moty, r_moty, s2p_isomorphism, m2r_isomorphism)
end
