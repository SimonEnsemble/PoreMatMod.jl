module MOFun

## packages that are used by this module
using PorousMaterials
using LightGraphs
using LinearAlgebra
using Printf
using DataFrames
using Moiety, Ullmann
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
function read_fragment_from_xyz(name::String)
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


@doc raw"""
Wraps find_subgraph_isomorphisms for convenient substructure searching
"""
function substructure_search(find_moiety::Crystal, parent_structure::Crystal)
	configurations = Ullmann.find_subgraph_isomorphisms(find_moiety.bonds,
											  find_moiety.atoms.species,
											  parent_structure.bonds,
											  parent_structure.atoms.species)
	return configurations
end



export
	# MOFfun.jl
	remove_extension, read_fragment_from_xyz, functionalize_mof,
	choose_side, substructure_search

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
