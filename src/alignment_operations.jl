####
# functions used for alignment of the Fragment to the AromaticRing
####

"""
Returns a formatted array with the ids of the atoms needed for the specific type of arene substtution.

Side two is functionalized by default.

ids = [atom_1, anchor, atom_2]

# Notation: 
- `anchor` referes to the atom on the aromatic ring that the functional group is going to be bonded to.
- `atom_1` is the atom on the aromatic ring bonded in the direction of the ipso atom.
- `atom_2` is the atom on the aromatic ring bonded in the direction of the para atom.
"""
function crystal_aro_triplet_ids(ring::AromaticRing, substitution_position::String, which_side::Int=2)
    if substitution_position == "meta"
        return [ring.sides[which_side].ortho.id_C, ring.sides[which_side].meta.id_C, ring.id_C_para]
    elseif substitution_position == "ortho"
        return [ring.id_ipso, ring.sides[which_side].ortho.id_C, ring.sides[which_side].meta.id_C]
    else
        error("substitution pattern not supported!")
    end
end

function fragment_aro_triplet_ids(fragment::Crystal)
    ids_C_aro = findall(fragment.atoms.species .== :C_aro)
    @assert length(ids_C_aro) == 2
    ids_C_aro_R = findall(fragment.atoms.species .== :C_aro_R)
    @assert length(ids_C_aro_R) == 1
    ids_aro_triplets = [ids_C_aro[1], ids_C_aro_R[1], ids_C_aro[2]]
    @assert fragment.atoms.species[ids_aro_triplets] == [:C_aro, :C_aro_R, :C_aro]
    return ids_aro_triplets
end

"""
Return the id of the atom bonded to :C_aro_R on the fragment.
"""
function fragment_aro_R_id(fragment::Crystal)
    frag_ids = 1:fragment.atoms.n
    id_frag_r = filter(id -> (have_neighbor(fragment, id, :C_aro_R) 
            && fragment.atoms.species[id] != :C_aro),
            frag_ids)
    @assert length(id_frag_r) == 1
    return id_frag_r[1]
end

"""
Find the geometric center of the atoms used for the substitution.
"""
function center_of_aro_triplet(crystal::Crystal, aro_triplet_ids::Array{Int, 1})
    @assert length(aro_triplet_ids) == 3
    return sum(crystal.atoms.coords.xf[:, aro_triplet_ids], dims=2)[:] / 3
end

"""
Return an array with the translated (if necessary), fractional coordinates of the `aro_triplet_crystal` 
to be used for alignment. We use the `anchor` atom as the reference atom for translation since it is the 
atom that will be functionalized. Translation is only necessary if the atoms in the aro_triplet_framewok
are split across periodic boundary conditions.
"""
function triplet_locality(crystal::Crystal, ids_aro_triplet::Array{Int64,1})
    # make sure atoms in triple are close to each other
    xf_aro_triplet = zeros(3, length(ids_aro_triplet))
    for (i, id_aro_triplet) in enumerate(ids_aro_triplet)
        # we want to use the anchor atoms as the reference position
        dxf = crystal.atoms.coords.xf[:, id_aro_triplet] - crystal.atoms.coords.xf[:, ids_aro_triplet[2]]
        nearest_image!(dxf)
        xf_aro_triplet[:, i] = dxf
    end
    return xf_aro_triplet
end

"""
Align fragment with the correct part of the aromatic ring using orthogonal procrustes
"""
function align_fragment(ring::AromaticRing, crystal::Crystal, 
        fragment::Crystal, substitution_position::String; which_side::Int=2)
    # find ids of [:C_aro, :C_aro_R, :C_aro]
    ids_aro_triplet_crystal = crystal_aro_triplet_ids(ring, substitution_position, which_side)
    ids_aro_triplet_fragment = fragment_aro_triplet_ids(fragment)
    
    # get the fractional coordinates of the atoms in the aro_triplet_crystal
    # and account for periodic boundary conditions
    xf_aro_triplet_crystal = triplet_locality(crystal, ids_aro_triplet_crystal)
    
    # get geometric center of "ghost" triplet in fractional coords
    xf_image_center = sum(xf_aro_triplet_crystal, dims=2)[:] / 3
    
    # determine the total coordinate tranlation of the system 
    #    to be added back in post-rotation
    # Fractional coords
    xf_offset = xf_image_center + crystal.atoms.coords.xf[:, ids_aro_triplet_crystal[2]] 
	 
    # find centers of aro-triplets for the fragment in fractional coords
    xf_center_fragment = center_of_aro_triplet(fragment, ids_aro_triplet_fragment)
    
    # translate centers of aro-triplets to the origin in fractional coords
    fragment.atoms.coords.xf .= fragment.atoms.coords.xf .- xf_center_fragment 
    xf_image_centered_crystal = xf_aro_triplet_crystal .- xf_image_center 
    
    # find rotation matrix to apply to the fragment to align its aro-triplet
    #   with the aro-triplet of the crystal.
    # set up orthogonal Procrustes
    A = fragment.box.f_to_c * fragment.atoms.coords.xf[:, ids_aro_triplet_fragment] # Cartesian
    B = crystal.box.f_to_c * xf_image_centered_crystal # Cartesian
    
    # find rotation matrix via orthogonal Procrustes in Cartesian space
    F = svd(A * B')
    R = F.V * F.U' # rotation matrix
    
    # apply transformation to the fragment
    #  rotate fragment about origin (it is currently centered)
    #  then translate it to the appropriate position in the crystal.
    aligned_coords = Atoms{Cart}(length(fragment.atoms.species),
        fragment.atoms.species, 
        Cart(R * fragment.box.f_to_c * fragment.atoms.coords.xf .+ crystal.box.f_to_c * xf_offset))
    
    # now convert the aligned atom positions into fractional coords w.r.t the crystal box
    aligned_coords = Frac(aligned_coords, crystal.box)    

    return Crystal(remove_extension(crystal) * "_aligned_" * fragment.name, 
        crystal.box,
        aligned_coords,
        Charges{Frac}(0),
        fragment.bonds,
        crystal.symmetry
    )
end
