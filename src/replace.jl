"""
    child = replace(parent, query => replacement)

Generates a `child` crystal structure by (i) searches the `parent` crystal structure for subgraphs that match the `query` then 
(ii) replaces the substructures of the `parent` matching the `query` fragment with the `replacement` fragment.

Equivalent to calling `substructure_replace(query ∈ parent, replacement)`.

Accepts the same keyword arguments as [`substructure_replace`](@ref).
"""
replace(p::Crystal, pair::Pair; kwargs...) = substructure_replace(pair[1] ∈ p, pair[2]; kwargs...)


"""
    alignment = Alignment(rot::Matrix{Float64}, shift_1::Vector{Float64}, shift_2::Vector{Float64}, err::Float64)

Data structure for tracking alignment in substructure find/replace operations.
"""
struct Alignment
	rot::Matrix{Float64}
	# before rotation
	shift_1::Vector{Float64}
	# after rotation
	shift_2::Vector{Float64}
	# error
	err::Float64
end


# Gets the rotation matrix for aligning the replacement moiety onto a subset (isomorphic to masked query) of parent atoms
function r2p_op(replacement::Crystal, parent::Crystal, r2p_isom::Dict{Int,Int})
    @assert replacement.atoms.n ≥ 3 && parent.atoms.n ≥ 3 "Parent and replacement must each be at least 3 atoms for SVD alignment."
    # get matrix A (replacement fragment coordinates)
    A = @views replacement.box.f_to_c * replacement.atoms[[r for (r,p) in r2p_isom]].coords.xf
    # prepare parent subset
    parent_subset = deepcopy(parent.atoms[[p for (r,p) in r2p_isom]].coords.xf)
    adjust_for_pb!(parent_subset)
    parent_subset_center = geometric_center(parent_subset)
    # get matrix B (parent subset coordinates)
    B = parent.box.f_to_c * (parent_subset .- parent_subset_center)
    # solve the SVD
    F = svd(A * B')
    # return rotation matrix
    return F.V * F.U'
end


# Transforms replacement according to two rotation matrices and a translational offset
function xform_replacement(replacement::Crystal, rot_r2p::Matrix{Float64}, parent_subset_center::Vector{Float64}, parent::Crystal)::Crystal
    # put replacement into cartesian space
    atoms = Atoms{Cart}(length(replacement.atoms.species), replacement.atoms.species, Cart(replacement.atoms.coords, replacement.box))
    # rotate replacement to align with parent_subset
    atoms.coords.x[:,:] = rot_r2p * atoms.coords.x
    # translate to align with original parent center
    atoms.coords.x .+= parent.box.f_to_c * parent_subset_center
    # cast atoms back to Frac
    xrm = Crystal(replacement.name, parent.box, Atoms{Frac}(length(atoms.species),
        atoms.species, Frac(atoms.coords, parent.box)), Charges{Frac}(0))
    # restore bonding network
    for e in edges(replacement.bonds)
        add_edge!(xrm.bonds, src(e), dst(e))
    end
    return xrm
end


# tracks which bonds need to be made between the parent and array
# of transformed replacement's (xrms) along with the new fragments
function accumulate_bonds!(bonds::Array{Tuple{Int,Int}}, q2p_isom::Dict{Int,Int},
    parent::Crystal, q_unmasked2r_isom::Union{Dict{Int,Int},Nothing}, xrm::Union{Crystal,Nothing}, count_xrms::Int)
    # skip bond accumulation for null replacement
    if q_unmasked2r_isom == Int[] || isnothing(q_unmasked2r_isom)
        return
    end
    # bonds between new fragments and parent
    # loop over q2p_isom
    for (s, (q, p)) in enumerate(q2p_isom)
        # in case the replacement moiety is smaller than the search moiety
        if s > length(q_unmasked2r_isom)
            break
        end
        # find neighbors of parent_subset atoms
        n = Graphs.neighbors(parent.bonds, p)
        # loop over neighbors
        for nᵢ in n
            # if neighbor not in q2p_isom, must bond it to replacement replacement
            # of parent_subset atom
            if !(nᵢ ∈ values(q2p_isom))
                # ID the atom in replacement
                r = q_unmasked2r_isom[s]
                # add the index offset
                r += parent.atoms.n + (count_xrms - 1) * xrm.atoms.n
                # push bond to array
                push!(bonds, (nᵢ, r))
            end
        end
    end
    # new fragment bonds
    # calculate indices in new xtal
    offset = (count_xrms - 1) * xrm.atoms.n + parent.atoms.n
    for e in edges(xrm.bonds) # loop over fragment edges
        push!(bonds, (offset + src(e), offset + dst(e)))
    end
end


# align the replacement onto the parent subgraph using the mappings between q, p, q′, and r
function align(q_unmasked2r_isom::Dict{Int,Int}, parent::Crystal, replacement::Crystal, r2p_isom::Dict{Int, Int}, q_unmasked2p_isom::Dict{Int,Int})
    # find parent subset
    parent_subset = deepcopy(parent[[q_unmasked2p_isom[q] for q in 1:length(q_unmasked2p_isom)]])
    # adjust coordinates for periodic boundaries
    adjust_for_pb!(parent_subset)
    # record the center of parent_subset so we can translate back later
    parent_subset_center = geometric_center(parent_subset)
    # shift all replacement nodes according to center of isomorphic subset
    replacement′ = deepcopy(replacement)
    replacement′.atoms.coords.xf .-= geometric_center(replacement[[q_unmasked2r_isom[q] for q in 1:length(q_unmasked2r_isom)]])
    write_cif(replacement′, "shifted_replacement.cif")
    # get OP rotation matrix to align replacement onto parent
    rot_r2p = r2p_op(replacement′, parent, r2p_isom)
    # transform replacement by rot_r2p, and parent_subset_center (this is now a potential crystal to add)
    xrm = xform_replacement(replacement′, rot_r2p, parent_subset_center, parent)
    # calculate the alignment error
    alignment_error = rmsd(xrm.atoms.coords.xf[:, [q_unmasked2r_isom[q] for q in 1:length(q_unmasked2r_isom)]], parent.atoms.coords.xf[:, [q_unmasked2p_isom[q] for q in 1:length(q_unmasked2p_isom)]])
    return xrm, alignment_error
end


# generates data for effecting a series of replacements
function build_replacement_data(configs::Vector{Tuple{Int,Int}}, q_in_p::Search, replacement::Crystal)
    parent = q_in_p.parent
    query = q_in_p.query
    # which atoms from query are in replacement?
    q_unmasked2q = [i for i in 1:query.atoms.n if !(i ∈ r_group_indices(query))]
    # BitArray for identifying atoms as masked (false) or unmasked (true)
    not_masked = map(species -> ! occursin(rc[:r_tag], String(species)), query.atoms.species) .== true
    # get isomrphism between query/mask and replacement
    q_unmasked_in_r = substructure_search(query[not_masked], replacement, assertion_override=true)
    # check for good input
    if nb_locations(q_unmasked_in_r) > 1
        @warn "There should be a single subset of the replacement isomorphic to the unmasked query!"
    end
    # take arbitrary isom, if any
    q_unmasked2r_isom = nb_isomorphisms(q_unmasked_in_r) ≠ 0 ? q_unmasked_in_r.isomorphisms[1][1] : nothing
    if !isnothing(q_unmasked2r_isom)
        q2r_isom′ = Dict([q => q_unmasked2r_isom[m] for (m, q) in enumerate(q_unmasked2q)])
    end
    # containers for optimal results
    xrms = Crystal[]
    del_atoms = Int[]
    bonds = Tuple{Int,Int}[] # tuple (i,j) encodes a bond in the child crystal
    # parent bonds
    for e in edges(parent.bonds) # loop over parent structure bonds
        push!(bonds, (src(e), dst(e))) # preserve them
    end
    # generate transformed replace moiety (xrm), ID which atoms to delete,
    # and accumulate list of bonds for each replacement configuration
    for config in configs
        loc = config[1]
        ori = config[2] == 0 ? [1:nb_ori_at_loc(q_in_p)[loc]...] : [config[2]]
        xrm = nothing
        alignment_err = Inf
        q2p_isom = nothing
        for ori′ in ori # choose best r2p by evaluating RMSD for all possibilities
            #### all variables w/ prime (e.g. ori′) are local to loop iteration
            # pull up specific isomorphism from query to parent
            q2p_isom′ = q_in_p.isomorphisms[loc][ori′]
            # orthog. Procrustes
            if !isnothing(q_unmasked2r_isom)
                # determine mapping r2p
                r2p_isom = Dict([r => q2p_isom′[q] for (q,r) in q2r_isom′])
                # determine isomorphism from masked query to parent
                q_unmasked2p_isom = Dict([q => q2p_isom′[q] for q in q_unmasked2q])
                # check error and keep xrm & q_unmasked2r_isom if better than previous best error
                xrm′, alignment_err′ = align(q_unmasked2r_isom, parent, replacement, r2p_isom, q_unmasked2p_isom)
                if alignment_err′ < alignment_err
                    alignment_err = alignment_err′
                    xrm = xrm′
                    q2p_isom = q2p_isom′
                end
            else
                q2p_isom = q2p_isom′ # the replacement is `nothing` (or invalid)
            end
        end
        # add optimal xrm to array
        if !isnothing(xrm)
            push!(xrms, xrm)
        end
        # push obsolete atoms to array
        for x in values(q2p_isom)
            push!(del_atoms, x) # this can create duplicates; remove them later
        end
        # clean up del_atoms
        del_atoms = unique(del_atoms)
        # accumulate bonds
        accumulate_bonds!(bonds, q2p_isom, parent, q_unmasked2r_isom, xrm, length(xrms))
    end
    return xrms, del_atoms, bonds
end


function aligned_replacement(replacement::Crystal, parent::Crystal, r2p_alignment::Alignment)
    # put replacement into cartesian space
    atoms_r = Cart(replacement.atoms, replacement.box)
    # rotate replacement to align with parent_subset
    atoms_r.coords.x[:, :] = r2p_alignment.rot * (atoms_r.coords.x .+ r2p_alignment.shift_1) .+ r2p_alignment.shift_2
    # cast atoms back to Frac
    return Crystal(replacement.name, parent.box, Frac(atoms_r, parent.box), Charges{Frac}(0))
end


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


# Internal method for performing substructure replacements
function __substructure_replace(q_in_p::Search, replacement::Crystal, configs::Array{Tuple{Int,Int}}, new_xtal_name::String)::Crystal

end

function _substructure_replace(q_in_p::Search, replacement::Crystal, configs::Array{Tuple{Int,Int}}, new_xtal_name::String)::Crystal ## old version (deleteme)
    parent = q_in_p.parent
    # configs must all be unique
    @assert length(configs) == length(unique(configs)) "configs must be unique"
    # mutation guard
    replacement = deepcopy(replacement)
    # if there are no replacements to be made, just return the parent
    if nb_isomorphisms(q_in_p) == 0
        @warn "No replacements to be made."
        return parent
    end
    # loop over configs to build replacement data
    xrms, del_atoms, bonds = build_replacement_data(configs, q_in_p, replacement)
    # append temporary crystals to parent
    atoms = xrms == Crystal[] ? parent.atoms : parent.atoms + sum([xrm.atoms for xrm ∈ xrms if xrm.atoms.n > 0])
    xtal = Crystal(new_xtal_name, parent.box, atoms, Charges{Frac}(0))
    # create bonds from tuple array
    for (i, j) ∈ bonds
        add_edge!(xtal.bonds, i, j)
    end
    # correct for periodic boundaries
    wrap!(xtal)
    # slice to final atoms/bonds
    new_xtal = xtal[[x for x ∈ 1:xtal.atoms.n if !(x ∈ del_atoms)]]
    # calculate bond attributes
    for bond ∈ edges(new_xtal.bonds)
        dist = distance(new_xtal.atoms, new_xtal.box, src(bond), dst(bond), true)
        cross_pb = dist != distance(new_xtal.atoms, new_xtal.box, src(bond), dst(bond), false)
        set_props!(new_xtal.bonds, bond, Dict(:distance => dist, :cross_boundary => cross_pb))
    end
    # handle symmetry
    if parent.symmetry.is_p1
        return new_xtal
    else
        @warn "Copying symmetry rules for non-P1 crystal"
        return Crystal(new_xtal.name, new_xtal.box, new_xtal.atoms, new_xtal.charges, new_xtal.bonds, parent.symmetry)
    end
end


@doc raw"""
    child = substructure_replace(search, replacement; random=false, nb_loc=0, loc=Int[], ori=Int[], name="new_xtal", verbose=false, remove_duplicates=false, periodic_boundaries=true)

Replace the substructures of `search.parent` matching the `search.query` fragment with the `replacement` fragment,
at locations and orientations specified by the keyword arguments `random`, `nb_loc`, `loc`, and `ori`.
Default behavior is to effect replacements at all "hit" locations in the parent structure and, at each location,
choose the orientation giving the optimal (lowest error) spatial aligment.
                
Returns a new `Crystal` with the specified modifications (returns `search.parent` if no replacements are made).

# Arguments
- `search::Search` the `Search` for a substructure moiety in the parent crystal
- `replacement::Crystal` the moiety to use for replacement of the searched substructure
- `random::Bool` set `true` to select random replacement orientations
- `nb_loc::Int` assign a value to select random replacement at `nb_loc` random locations
- `loc::Array{Int}` assign value(s) to select specific locations for replacement.  If `ori` is not specified, replacement orientation is random.
- `ori::Array{Int}` assign value(s) when `loc` is assigned to specify exact configurations for replacement. `0` values mean the configuration at that location should be selected for optimal alignment with the parent.
- `name::String` assign to give the generated `Crystal` a name ("new_xtal" by default)
- `verbose::Bool` set `true` to print console messages about the replacement(s) being performed
- `remove_duplicates::Bool` set `true` to automatically combine overlapping atoms of the same species in generated structures.
- `reinfer_bonds::Bool` set `true` to re-infer bonds after producing a structure
- `periodic_boundaries::Bool` set `false` to disable periodic boundary conditions when checking for atom duplication or re-inferring bonds
"""
function substructure_replace(search::Search, replacement::Crystal; random::Bool=false,
    nb_loc::Int=0, loc::Array{Int}=Int[], ori::Array{Int}=Int[], name::String="new_xtal", verbose::Bool=false,
    remove_duplicates::Bool=false, periodic_boundaries::Bool=true, reinfer_bonds::Bool=false)::Crystal
    # replacement at all locations (default)
    if nb_loc == 0 && loc == Int[] && ori == Int[]
        nb_loc = nb_locations(search)
        loc = [1:nb_loc...]
        if random
            ori = [rand(1:nb_ori_at_loc(search)[i]) for i in loc]
            if verbose
                @info "Replacing" q_in_p=search r=replacement.name mode="random ori @ all loc"
            end
        else
            ori = zeros(Int, nb_loc)
            if verbose
                @info "Replacing" q_in_p=search r=replacement.name mode="optimal ori @ all loc"
            end
        end
    # replacement at nb_loc random locations
    elseif nb_loc > 0 && ori == Int[] && loc == Int[]
        loc = sample([1:nb_locations(search)...], nb_loc, replace=false)
        if random
            ori = [rand(1:nb_ori_at_loc(search)[i]) for i in loc]
            if verbose
                @info "Replacing" q_in_p=search r=replacement.name mode="random ori @ $nb_loc loc"
            end
        else
            ori = zeros(Int, nb_loc)
            if verbose
                @info "Replacing" q_in_p=search r=replacement.name mode="optimal ori @ $nb_loc loc"
            end
        end
    # specific replacements
    elseif ori ≠ Int[] && loc ≠ Int[]
        @assert length(loc) == length(ori) "one orientation per location"
        nb_loc = length(ori)
        if verbose
            @info "Replacing" q_in_p=search r=replacement.name mode="loc: $loc\tori: $ori"
        end
    # replacement at specific locations
    elseif loc ≠ Int[]
        nb_loc = length(loc)
        if random
            ori = [rand(1:nb_ori_at_loc(search)[i]) for i in loc]
            if verbose
                @info "Replacing" q_in_p=search r=replacement.name mode="random ori @ loc: $loc"
            end
        else
            ori = zeros(Int, nb_loc)
            if verbose
                @info "Replacing" q_in_p=search r=replacement.name mode="optimal ori @ loc: $loc"
            end
        end
    end

    # generate configuration tuples (location, orientation)
    configs = Tuple{Int,Int}[(loc[i], ori[i]) for i in 1:nb_loc]
    # process replacements
    child = _substructure_replace(search, replacement, configs, name)

    if remove_duplicates
        child = Crystal(child.name, child.box, 
            Xtals.remove_duplicates(child.atoms, child.box, periodic_boundaries),
            Xtals.remove_duplicates(child.charges, child.box, periodic_boundaries)
        )
    end

    if reinfer_bonds
        remove_bonds!(child)
        infer_bonds!(child, periodic_boundaries)
    end

    return child
end

substructure_replace(search::Search, replacement::Nothing; kwargs...) = substructure_replace(search, moiety(nothing); kwargs...)
