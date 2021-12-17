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


struct Installation
    aligned_replacement::Crystal
    q2p::Dict{Int, Int}
    r2p::Dict{Int, Int}
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
    parent_substructure = deepcopy(parent[[p for (r, p) in r2p]])
    conglomerate!(parent_substructure)
    atoms_p = Cart(parent_substructure.atoms, parent_substructure.box)
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


function conglomerate!(parent_substructure::Crystal)
    # snip the cross-PB bonds
    bonds = deepcopy(parent_substructure.bonds)
    if length(connected_components(bonds)) > 1
        @warn "# connected components in parent substructure > 1. assuming the substructure does not cross the periodic boundary..."
        return
    end
    drop_cross_pb_bonds!(bonds)

    # find connected components of bonding graph without cross-PB bonds
    #    these are the components split across the boundary
    conn_comps = connected_components(bonds)
    
    # if substructure is entireline in the unit cell, it's already conglomerated :)
    if length(conn_comps) == 1
        return
    end

    # we wish to shift all connected components to a reference component,
    #   defined to be the largest component for speed.
    conn_comps_shifted = [false for c = 1:length(conn_comps)]
    ref_comp_id = argmax(length.(conn_comps))
    conn_comps_shifted[ref_comp_id] = true  # consider it shifted.
    
    # has atom p been shifted?
    function shifted_atom(p::Int)
        # loop over all connected components that have been shifted
        for conn_comp in conn_comps[conn_comps_shifted]
            # if parent substructure atom in this, yes!
            if p in conn_comp
                return true
            end
        end
        # reached this far, atom p is not in component that has been shifted.
        return false
    end
    
    # to which component does atom p belong?
    function find_component(p::Int)
        for c = 1:length(conn_comps)
            if p in conn_comps[c]
                return c
            end
        end
    end
    
    # until all components have been shifted to the reference component...
    while ! all(conn_comps_shifted)
        # loop over cross-PB edges in the parent substructure
        for ed in edges(parent_substructure.bonds)
            if get_prop(parent_substructure.bonds, ed, :cross_boundary)
                # if one edge belongs to unshifted component and another belogs to any component that has been shifted...
                if     shifted_atom(ed.src) && ! shifted_atom(ed.dst)
                    p_ref, p = ed.src, ed.dst
                elseif shifted_atom(ed.dst) && ! shifted_atom(ed.src)
                    p_ref, p = ed.dst, ed.src
                else
                    continue # both are shifted or both are unshifted. ignore this cross-PB edge
                end
                # here's the unshifted component we will shift next, to be next to the shifted components.
                comp_id = find_component(p)
                # find displacement vector for this cross-PB edge.
                dx = parent_substructure.atoms.coords.xf[:, p_ref] - parent_substructure.atoms.coords.xf[:, p]
                # get distance to nearest image
                n_dx = copy(dx)
                nearest_image!(n_dx)
                # shift all atoms in this component by this vector.
                for atom_idx in conn_comps[comp_id]
                    parent_substructure.atoms.coords.xf[:, atom_idx] .+= dx - n_dx
                end
                # mark that we've shifted this component.
                conn_comps_shifted[comp_id] = true
            end
        end
    end
    return
end


function aligned_replacement(replacement::Crystal, parent::Crystal, r2p_alignment::Alignment)
    # put replacement into cartesian space
    atoms_r = Cart(replacement.atoms, replacement.box)
    # rotate replacement to align with parent_subset
    atoms_r.coords.x[:, :] = r2p_alignment.rot * (atoms_r.coords.x .+ r2p_alignment.shift_1) .+ r2p_alignment.shift_2
    # cast atoms back to Frac
    return Crystal(replacement.name, parent.box, Frac(atoms_r, parent.box), Charges{Frac}(0), replacement.bonds, replacement.symmetry)
end


function effect_replacements(search::Search, replacement::Crystal, configs::Vector{Tuple{Int, Int}}, name::String)::Crystal
    nb_not_masked = sum(.! occursin.(rc[:r_tag], String.(search.query.atoms.species)))
    if replacement.atoms.n > 0
        q_unmasked_in_r = substructure_search(search.query[1:nb_not_masked], replacement)
        q2r = Dict([q => q_unmasked_in_r.isomorphisms[1][1][q] for q in 1:nb_not_masked])
    else
        q2r = Dict{Int, Int}()
    end
    
    installations = [optimal_replacement(search, replacement, q2r, loc_id, [ori_id]) for (loc_id, ori_id) in configs]
    
    child = install_replacements(search.parent, installations, name)

    # handle `missing` values in edge :cross_boundary attribute
    for edge in edges(child.bonds) # loop over edges
        # check if cross-boundary info is missing
        if ismissing(get_prop(child.bonds, edge, :cross_boundary))
            # check if bond crosses boundary
            distance_e = get_prop(child.bonds, edge, :distance) # distance in edge property
            dxa = Cart(Frac(child.atoms.coords.xf[:, src(edge)] - child.atoms.coords.xf[:, dst(edge)]), child.box) # Cartesian displacement
            distance_a = norm(dxa.x) # current euclidean distance by atom coords
            set_prop!(child.bonds, edge, :cross_boundary, !isapprox(distance_e, distance_a, atol=0.1))
        end
    end

    return child
end


function install_replacements(parent::Crystal, replacements::Vector{Installation}, name::String)::Crystal
    # create child w/o symmetry rules for sake of crystal addition
    child = Crystal(name, parent.box, parent.atoms, parent.charges, parent.bonds, Xtals.SymmetryInfo())

    obsolete_atoms = Int[] # to delete at the end

    # loop over replacements to install
    for installation in replacements
        replacement, q2p, r2p = installation.aligned_replacement, installation.q2p, installation.r2p
        #add into parent
        if replacement.atoms.n > 0
            child = +(child, replacement, check_overlap=false)
        end
        
        # reconstruct bonds
        for (r, p) in r2p # p is in parent_subst
            p_nbrs = neighbors(parent.bonds, p)
            for p_nbr in p_nbrs
                if ! (p_nbr in values(q2p)) # p_nbr not in parent_subst
                    # need bond nbr => r in child, where r is in replacement
                    e = (p_nbr, child.atoms.n - replacement.atoms.n + r)
                    # create edge
                    add_edge!(child.bonds, e)
                    # copy edge attributes from parent (:cross_boundary will need to be reassessed later)
                    set_props!(child.bonds, e[1], e[2], props(parent.bonds, p, p_nbr))
                    set_prop!(child.bonds, e[1], e[2], :cross_boundary, missing)
                end
            end
        end
        
        # accumulate atoms to delete
        obsolete_atoms = vcat(obsolete_atoms, values(q2p)...)
    end

    # delete obsolete atoms
    obsolete_atoms = unique(obsolete_atoms)
    keep_atoms = [p for p = 1:child.atoms.n if ! (p in obsolete_atoms)]
    child = child[keep_atoms]

    # restore symmetry rules
    child = Crystal(name, child.box, child.atoms, child.charges, child.bonds, parent.symmetry)

    # return result
    return child
end


function optimal_replacement(search::Search, replacement::Crystal, q2r::Dict{Int,Int}, loc_id::Int, ori_ids::Vector{Int})
    # unpack search arg
    isomorphisms, parent = search.isomorphisms, search.parent

    if q2r == Dict{Int, Int}() # "replace-with-nothing" operation
        q2p = isomorphisms[loc_id][1]
        r2p = Dict([0 => p for p in values(q2p)])
        return Installation(replacement, q2p, r2p)
    end

    if ori_ids == [0]
        ori_ids = [1:nb_ori_at_loc(search)[loc_id]...]
    end

    # loop over ori_ids to find best r2p_alignment
    r2p_alignment = Alignment(zeros(1,1), [0.], [0.], Inf)
    best_ori = 0
    best_r2p = Dict{Int, Int}()
    for ori_id in ori_ids
        # find r2p isom
        q2p = isomorphisms[loc_id][ori_id]
        r2p = Dict([r => q2p[q] for (q, r) in q2r])
        # calculate alignment
        test_alignment = get_r2p_alignment(replacement, parent, r2p)
        # keep best alignment and generating ori_id
        if test_alignment.err < r2p_alignment.err
            r2p_alignment = test_alignment
            best_ori = ori_id
            best_r2p = r2p
        end
    end

    opt_aligned_replacement = aligned_replacement(replacement, parent, r2p_alignment)
    
    # return the replacement modified according to r2p_alignment
    @assert ne(opt_aligned_replacement.bonds) == ne(replacement.bonds)
    return Installation(opt_aligned_replacement, isomorphisms[loc_id][best_ori], best_r2p)
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
    remove_duplicates::Bool=false, periodic_boundaries::Bool=true, reinfer_bonds::Bool=false, wrap::Bool=true)::Crystal
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

    # remove charges from parent
    if search.parent.charges.n > 0
        @warn "Dropping charges from parent."
        p = search.parent
        search = Search(
            Crystal(p.name, p.box, p.atoms, Charges{Frac}(0), p.bonds, p.symmetry),
            search.query,
            search.isomorphisms
        )
    end

    # generate configuration tuples (location, orientation)
    configs = Tuple{Int,Int}[(loc[i], ori[i]) for i in 1:nb_loc]
    # process replacements
    child = effect_replacements(search, replacement, configs, name)

    if remove_duplicates
        child = Crystal(child.name, child.box, 
            Xtals.remove_duplicates(child.atoms, child.box, periodic_boundaries),
            Xtals.remove_duplicates(child.charges, child.box, periodic_boundaries)
        )
    end

    if wrap
        # wrap coordinates
        wrap!(child.atoms.coords)
        # check :cross_boundary edge attributes
        for edge in edges(child.bonds) # loop over edges
            distance_e = get_prop(child.bonds, edge, :distance) # distance in edge property
            dxa = Cart(Frac(child.atoms.coords.xf[:, src(edge)] - child.atoms.coords.xf[:, dst(edge)]), child.box) # Cartesian displacement
            distance_a = norm(dxa.x) # current euclidean distance by atom coords
            set_prop!(child.bonds, edge, :cross_boundary, !isapprox(distance_e, distance_a, atol=0.1))
        end
    end

    if reinfer_bonds
        remove_bonds!(child)
        infer_bonds!(child, periodic_boundaries)
    end

    return child
end

substructure_replace(search::Search, replacement::Nothing; kwargs...) = substructure_replace(search, moiety(nothing); kwargs...)
