
## Globals
# set to true for debugging output of xyz files
DEBUG = false

## Structs, conversions, and helpers
import Base.convert
# for encoding the location and orientation indices within a group of search results
struct Configuration
    location::Int
    orientation::Int
end
# hack to fix weird bug
convert(Configuration, (i,j)) = Configuration(i,j)
# for describing the search query ("find x in y")
struct SubstructureQuery
    substructure::Crystal
    parent_structure::Crystal
end
# hack for convenience
convert(SubstructureQuery,
	(a,b)::Tuple{Crystal, Crystal}) = SubstructureQuery(a,b)
# hack to display queries nicely
function Base.show(io::IO, sq::SubstructureQuery)
	print("$(sq.substructure.name) ∈ $(sq.parent_structure.name)")
end
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
## Helper functions
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
import Base.∈
(∈)(s::Crystal, g::Crystal) = substructure_search(s, g)

# Helper for making debugging .xyz's
import PorousMaterials.write_xyz
write_xyz(xtal::Crystal, name::String) = write_xyz(Cart(xtal.atoms, xtal.box), name)

# writes a debuggin output of the input Atoms as xyz
function debugxyz(atoms::Atoms{Cart}, str::String)
	if DEBUG
		@debug "debugxyz" str
		try
			write_xyz(atoms, "output/debug/" * str)
		catch
			mkdir("output/debug")
			write_xyz(atoms, "output/debug/" * str)
		end
	end
end
debugxyz(xtal::Crystal, str::String) = debugxyz(Cart(xtal.atoms, xtal.box), str)
debugxyz(t::Tuple{Crystal, String}) = debugxyz(t[1], t[2])
debugxyz(a::Array{Tuple{Crystal, String}}) = [debugxyz(x[1], x[2]) for x in a]

@doc raw"""
Searches for a substructure within a `Crystal` and returns a `SubstructureSearchResults`
struct containing all identified subgraph isomorphisms.
"""
function substructure_search(s_moty::Crystal, xtal::Crystal)::SubstructureSearchResults
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
	num_orientations = [length(results[location_set]) for location_set in keys(results)]
	results = [result for result in values(results)]
	# repack as array of structs
	result_struct_array = Array{SubstructureSearchResult}([])
    for (location, isomorphisms) ∈ enumerate(results)
        for (orientation, isomorphism) ∈ enumerate(results[location])
            push!(result_struct_array, SubstructureSearchResult((s_moty, xtal),
				(location, orientation), isomorphism))
        end
    end
    return SubstructureSearchResults(result_struct_array, length(result_struct_array),
		length(keys(results)), num_orientations)
end

# Translates all atoms in xtal such that Atoms[1] is @ origin, w/ correction for
# the Crystal's periodic boundaries.
function adjust_for_pb!(xtal::Crystal)
	origin_offset = deepcopy(xtal.atoms.coords.xf[:, 1])
    for i in 1:xtal.atoms.n
		 # rebuild the image (fixing first node @ origin)
		dxf = xtal.atoms.coords.xf[:, i] .- origin_offset
		# nearest_image! expects points to be within same or adjacent unit cells
		@assert all(abs.(dxf) .< 2) "Invalid xf coords in $(xtal.name)"
        nearest_image!(dxf)
        xtal.atoms.coords.xf[:, i] .= dxf
    end
end

# Performs orthogonal Procrustes on correlated point clouds A and B
function orthogonal_procrustes(A::Array{Float64,2}, B::Array{Float64,2})::Array{Float64,2}
    F = svd(A * B')
    return F.V * F.U'
end

# Gets the s_moty-to-xtal rotation matrix
function s2p_op(s_moty::Crystal, xtal::Crystal,	s2p_isomorphism::Array{Int})::Array{Float64}
	A = s_moty.box.f_to_c * s_moty.atoms.coords.xf
	B = xtal.box.f_to_c * xtal.atoms.coords.xf[:, s2p_isomorphism]
	return orthogonal_procrustes(A, B)
end

# Gets the r_moty-to-s_mask rotation matrix
function r2m_op(r_moty::Crystal, s_moty::Crystal, m2r_isomorphism::Array{Int},
		s_mask_atoms::Atoms)::Array{Float64}
	A = r_moty.box.f_to_c * r_moty.atoms[m2r_isomorphism].coords.xf
	B = s_moty.box.f_to_c * s_mask_atoms.coords.xf
	return orthogonal_procrustes(A, B)
end

# Transforms r_moty according to two rotation matrices and a translational offset
function xform_r_moty(r_moty::Crystal, rot_r2m::Array{Float64,2}, rot_s2p::Array{Float64,2},
		xtal_offset::Array{Float64}, xtal::Crystal)::Crystal
	# put r_moty into cartesian space
	atoms = Atoms{Cart}(length(r_moty.atoms.species), r_moty.atoms.species,
		Cart(r_moty.atoms.coords, r_moty.box))
	# transformation 1: rotate r_moty to align with s_moty
	atoms.coords.x[:,:] = rot_r2m * atoms.coords.x
	debugxyz(atoms, "12_rotation1")
	# transformation 2: rotate to align with xtal_subset
	atoms.coords.x[:,:] = rot_s2p * atoms.coords.x
	debugxyz(atoms, "13_rotation2")
	# transformation 3: translate to align with original xtal center
	atoms.coords.x .+= xtal.box.f_to_c * xtal_offset
	debugxyz(atoms, "14_translation2")
	# cast atoms back to Frac
	return Crystal(r_moty.name, xtal.box, Atoms{Frac}(length(atoms.species), atoms.species, Frac(atoms.coords, xtal.box)), Charges{Frac}(0))
end

# shifts coordinates to make the geometric center of the point cloud coincident w/ the origin
function center_on_self!(xtal::Crystal)
	xtal.atoms.coords.xf .-= geometric_center(xtal)
end

# returns a Crystal with a specific replacement made
function effect_replacement(xtal::Crystal, s_moty::Crystal,	r_moty::Crystal,
		s2p_isomorphism::Array{Int}, m2r_isomorphism)::Crystal
	# prevent input mutation and slice off the target subset of the parent structure
	s_moty, r_moty, xtal_subset = deepcopy.([s_moty, r_moty, xtal[s2p_isomorphism]])
	debugxyz.([(xtal, "01_xtal"), (r_moty, "02_r_moty"), (s_moty, "03_s_moty")])
	xtal_subset_center = geometric_center(xtal_subset) # where to return post-alignment
	# adjust the coordinates for periodic boundaries
    adjust_for_pb!(xtal_subset)
	# shift to align centers at origin
	center_on_self!.([xtal_subset, s_moty])
	debugxyz.([(xtal_subset, "04_adjusted_xtal_subset"), (s_moty, "05_adjusted_s_moty")])
    # determine s_mask (which atoms from s_moty are NOT in r_moty?)
    s_mask_atoms = s_moty.atoms[[i for i in 1:s_moty.atoms.n if !(i ∈ r_group_indices(s_moty))]]
	debugxyz(Cart(s_mask_atoms, s_moty.box), "08_search_moiety_mask")
    # shift all r_moty nodes according to center of isomorphic subset
    r_moty.atoms.coords.xf .-= geometric_center(r_moty[m2r_isomorphism])
	debugxyz(r_moty, "11_shifted_r_moty")
    # do orthogonal Procrustes for s_moty-to-parent and mask-to-replacement alignments
	rot_s2p = s2p_op(s_moty, xtal, s2p_isomorphism)
	rot_r2m = r2m_op(r_moty, s_moty, m2r_isomorphism, s_mask_atoms)
    # transform r_moty according to rot_r2m, rot_s2p, and xtal_subset_center, align to box
    r_moty = xform_r_moty(r_moty, rot_r2m, rot_s2p, xtal_subset_center, xtal)
    # subtract s_moty isomorphic subset from xtal, add back r_moty
    xtal = Crystal(remove_extension(xtal) * "_find_" *
		MOFun.remove_path_prefix(s_moty.name) *	"_replace_" * r_moty.name, xtal.box,
		Frac(Cart(xtal.atoms[
			[i for i in 1:length(xtal.atoms.species) if !(i ∈ s2p_isomorphism)]],
			xtal.box) + Cart(r_moty.atoms, r_moty.box),	xtal.box), Charges{Frac}(0))
	debugxyz(xtal, "15_crystal")
    wrap!(xtal)
    infer_bonds!(xtal, true) # this needs to be changed to intelligent bond creation!
    return xtal
end

@doc raw"""
Finds the search moiety in the parent structure and replaces it at a location.
"""
function find_replace(s_moty::Crystal, r_moty::Crystal, xtal::Crystal;
		c::Int=1, s::SubstructureSearchResults=nothing)::Crystal
	if s == nothing
		s = s_moty ∈ xtal
	end
	mask = subtract_r_group(s_moty)
	infer_bonds!(mask, false)
	return effect_replacement(xtal, s_moty, r_moty,
		(s_moty ∈ xtal).results[c].isomorphism, (mask ∈ r_moty).results[1].isomorphism)
end
