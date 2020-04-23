####
# Structs used for aromatic rings
####
# these represent Carbon-Hydrogen pairs in the ring  
mutable struct CH
    id_C::Int
    id_H::Int
end

# on each side of the ring, there is an ortho and a meta position\
#  that is occupied by a carbon bonded to a Hydrogen
mutable struct RingSide
    ortho::CH
    meta::CH
end

# define the structure of the aromatic rings
mutable struct AromaticRing
    id_ipso::Int
    sides::Array{RingSide, 1}
    id_C_para::Int
end

# define an empty set for each kind of struct, so that we have a blank template that we can
# populate while building and labeling rings
empty_CH() = CH(-1, -1)
empty_ring_side() = RingSide(empty_CH(), empty_CH())
empty_ring() = AromaticRing(-1, [empty_ring_side(), empty_ring_side()], -1) 

# initialize empty ring struct (goes in main function?)
empty_ring()

####
# helper functions for filter
####
# given a cycle, check if it is an aromatic ring
# an aromatic ring is defined as a cycle whose atoms are all :C or `ipso_species`
function is_aromatic(crystal::Crystal, cyc::Array{Int, 1}, ipso_species::Symbol)
    for a in cyc
        if ! (crystal.atoms.species[a] in [:C, ipso_species])
            return false
        end
    end
    # if made it this far, only ipso and C in ring
    return true
end

## find unique aromatic cycles in the crystal
function find_aromatic_cycles(crystal::Crystal, ipso_species::Symbol, n::Int=6)
    # find all cycles less than or equal to length n
    all_cycles = simplecycles_limited_length(crystal.bonds, n, 10^6)
    
    # find all cycles that are aromatic and of length equal to n
    aro_cycles = filter(cyc -> (length(cyc) == n && is_aromatic(crystal, cyc, ipso_species)), all_cycles)
    
    # find unique cycles (above gives cycles in each direction, double counting)
    unique_aro_cycles = [] # let us fill this stack of cycles.
    for cyc in aro_cycles
        # flag to determine if the cycle is already in stack
        cyc_already_in_stack = false
        # if this cycle is already in our list of unique cycles, don't push it to the stack!
        for unique_cyc in unique_aro_cycles
            if issetequal(cyc, unique_cyc)
               cyc_already_in_stack = true
               break
            end
        end
        
        # if cycle wasn't in stack, push it to stack.
        if ! cyc_already_in_stack
            push!(unique_aro_cycles, cyc)
        end
    end
    
    return unique_aro_cycles
end

# is atom `i` in `crystal` `species`?
is_species(crystal::Crystal, i::Int, species::Symbol) = crystal.atoms.species[i] == species

## is atom `i` in `crystal` bonded to a `species`? 
function have_neighbor(crystal::Crystal, i::Int, species::Symbol)
    ids_neighbors = neighbors(crystal.bonds, i)
    return species in crystal.atoms.species[ids_neighbors]
end

# find the index of the ipso_species
function find_ipso_id(crystal::Crystal, r_species::Symbol, ipso_species::Symbol, cyc::Array{Int64, 1})
    # loop through atoms in the cycle tolook for the ipso atom
    for atom_id in cyc
        # is this atom the ipso species?
        if ! is_species(crystal, atom_id, ipso_species)
            continue
        end
        
        # find ids of neighbors, *excluding those already in the cycle*
        ids_neighbors = setdiff(neighbors(crystal.bonds, atom_id), cyc)
        
        # is this atom bonded to an R-species that is NOT in the cycle? then it is an ipso!
        for id_neighbor in ids_neighbors
            if is_species(crystal, id_neighbor, r_species)
                return atom_id
            end
        end
        # note: the ipso atom is not necessarily a unique choice, e.g. the C in IRMOF-1: there are two!
    end
    # if made it this far, then no ipso species was found
    error("no atom in the cycle was identified as an ipso species")
end

####
# filter function for constucting aromatic rings from cycles in crystal
####
function cycle_to_ring(crystal::Crystal, r_species::Symbol, ipso_species::Symbol, cyc::Array{Int64, 1})
    # find the ipso atom in the cycle
    id_ipso = find_ipso_id(crystal, r_species, ipso_species, cyc)
    
    # construct a ring that we'll populate.
    ring = AromaticRing(id_ipso, [empty_ring_side(), empty_ring_side()], -1)
    
    ###
    #   identify ortho carbons
    #      - must be a carbon atom
    #      - ipso atom is its neighbor
    #      - must be in cyc
    ###
    ids_C_ortho = filter(
        id -> is_species(crystal, id, :C) && 
              (ring.id_ipso in neighbors(crystal.bonds, id)), 
        cyc)
    @assert length(ids_C_ortho) == 2 "should be two ortho C's!"
    ring.sides[1].ortho.id_C, ring.sides[2].ortho.id_C = ids_C_ortho

    all_ids = [i for i = 1:crystal.atoms.n] # all atom ids
    for side in ring.sides
        ###
        #    identify H atoms bound to the ortho carbons
        #      - must be a H atom
        #      - must be bonded to the ortho C
        ###
        ids_H_ortho = filter(id -> is_species(crystal, id, :H) &&
                                   id in neighbors(crystal.bonds, side.ortho.id_C), 
                             all_ids)
        @assert length(ids_H_ortho) == 1 "should be only one H connected to C_ortho!"
        side.ortho.id_H = ids_H_ortho[1]

        ###
        #   identify meta carbons
        #      - must be a carbon atom
        #      - ortho atom is its neighbor
        #      - not the ipso atom
        #      - must be in cyc
        ###
        ids_C_meta = filter(id -> is_species(crystal, id, :C) && 
                                  id in neighbors(crystal.bonds, side.ortho.id_C) && 
                                  ring.id_ipso != id, 
                            cyc)
        @assert length(ids_C_meta) == 1 "should be only one C meta atom"
        side.meta.id_C = ids_C_meta[1]

        ###
        #    identify H atoms bound to the meta carbons
        #      - must be a H atom
        #      - must be bonded to the meta C
        ###
        ids_H_meta = filter(id -> is_species(crystal, id, :H) && 
                                  id in neighbors(crystal.bonds, side.meta.id_C), 
                            all_ids)
        @assert length(ids_H_meta) == 1 "should be only one H meta atom"
        side.meta.id_H = ids_H_meta[1]
    end
    
    ### 
    #   identify the para carbon
    #      - must be a carbon atom
    #      - must be bonded to *both* meta carbons
    #      - must be in cyc
    ###
    id_C_para = filter(id -> id in neighbors(crystal.bonds, ring.sides[1].meta.id_C) && 
                             id in neighbors(crystal.bonds, ring.sides[2].meta.id_C) && 
                             is_species(crystal, id, :C),
                       cyc)
    @assert length(id_C_para) == 1 "should be only one C para atom!"
    ring.id_C_para = id_C_para[1]
    
    ###
    # run sanity checks (tons of assert statements).
    # for example, the number of ortho carbons should be twice as large as the number of ipso atoms.
    ###
    @assert crystal.atoms.species[ring.id_C_para] == :C
    for i = 1:2
        @assert crystal.atoms.species[ring.sides[i].meta.id_C] == :C
        @assert crystal.atoms.species[ring.sides[i].meta.id_H] == :H
        @assert crystal.atoms.species[ring.sides[i].ortho.id_C] == :C
        @assert crystal.atoms.species[ring.sides[i].ortho.id_H] == :H
    end
    for side in ring.sides
        # ortho C connected to meta C
        @assert side.meta.id_C in neighbors(crystal.bonds, side.ortho.id_C)
        # H connected to C
        @assert side.meta.id_H in neighbors(crystal.bonds, side.meta.id_C)
        @assert side.ortho.id_H in neighbors(crystal.bonds, side.ortho.id_C)
    end
    
    return ring
end

"""
given the R-atom and ipso-atom species, identify aromatic rings in the structure and populate the ipso field.
returns an array of AromaticRings. 

For highly symetric aromatic rings (like those found in IRMOF-1) there is an ambiguity about which atom 
will be labeled `ipso`. By default it is the first atom in the cycle that satissfies the requirments to be ipso.
"""
function find_rings(crystal::Crystal, r_species::Symbol, ipso_species::Symbol; n::Int=6)
    # find all aromatic cyles
    cycles = find_aromatic_cycles(crystal, ipso_species, n)
    # loop through cycles, construct rings
    rings = AromaticRing[]
    for cyc in cycles
        ring = cycle_to_ring(crystal, r_species, ipso_species, cyc)
        push!(rings, ring)
    end
    return rings
end
