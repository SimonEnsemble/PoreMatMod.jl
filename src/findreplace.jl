## Structs
"""
    SearchTerms(parent, query)

Stores the `Crystal` inputs used to generate a `Search`
"""
struct SearchTerms
    parent::Crystal
    query::Crystal
end

Base.show(io::IO, q::SearchTerms) = print(io, q.query.name, " ∈ ", q.parent.name)


"""
    Search(search_terms, results)

Stores the `SearchTerms` used for a substructure search, and the results `DataFrame`
returned by carrying out the search.  Results are grouped by location in the
parent `Crystal` and can be examined using `nb_isomorphisms`, `nb_locations`,
and `nb_configs_at_loc`.  Subgraph isomorphisms are encoded like

    `isom = [7, 21, 9]`

where `isom[i]` is the index of the atom in `search.search.parent` corresponding
to atom `i` in `search.search.query` for the location and orientation specific
to `isom`.
"""
struct Search
    search::SearchTerms # the search terms (query ∈ parent)
    results
end

Base.show(io::IO, s::Search) = begin
    println(io, s.search)
    print(io, nb_isomorphisms(s), " hits in ", nb_locations(s), " locations.")
end


## Helpers
@doc raw"""
    nb_isomorphisms(search::Search)

Returns the number of isomorphisms found in the specified `Search`

# Arguments
- `search::Search` a substructure `Search` object
"""
function nb_isomorphisms(search::Search)::Int
    return sum([size(grp, 1) for grp in search.results])
end


@doc raw"""
    nb_locations(search::Search)

Returns the number of unique locations (collections of atoms) at which the
specified `Search` results contain isomorphisms.

# Arguments
- `search::Search` a substructure `Search` object
"""
function nb_locations(search::Search)::Int
    return length([size(grp, 1) for grp in search.results])
end


@doc raw"""
    nb_configs_at_loc(search)

Returns a array containing the number of isomorphic configurations at a given
location (collection of atoms) for which the specified `Search` results
contain isomorphisms.

# Arguments
- `search::Search` a substructure `Search` object
"""
function nb_configs_at_loc(search::Search)::Array{Int}
    return [size(grp, 1) for grp in search.results]
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
#    s ∈ g                   →    find the moiety in the crystal
#    [s1, s2] .∈ [g]         →    find each moiety in a crystal
#    s .∈ [g1, g2]           →    find the moiety in each crystal
#    [s1, s2] .∈ [g1, g2]    →    find each moiety in each crystal
(∈)(s::Crystal, g::Crystal) = substructure_search(s, g)


# extension of Base.replace to allow for simple and syntactically obvious find-replace
replace(p::Crystal, pair::Pair; kwargs...) = substructure_replace(pair[1] ∈ p, pair[2]; kwargs...)


# Helper for making .xyz's
write_xyz(xtal::Crystal, name::String) = Xtals.write_xyz(Cart(xtal.atoms, xtal.box), name)


"""
iso_structs = isomorphic_substructures(s::Search)::Crystal

Returns a crystal consisting of the atoms involved in subgraph isomorphisms in the search `s`
"""
function isomorphic_substructures(s::Search)::Crystal
    p = s.search.parent
    n = nb_locations(s)
    r = s.results
    is = sum([p[r[i].isomorphism[1]] for i in 1:n])
end


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


# Gets the query-to-parent rotation matrix
function s2p_op(query::Crystal, parent::Crystal)::Array{Float64,2}
    # query in Cartesian
    A = query.box.f_to_c * query.atoms.coords.xf
    # parent subset in Cartesian
    B = parent.box.f_to_c * parent.atoms.coords.xf
    # get rotation matrix
    return orthogonal_procrustes(A, B)
end


# Gets the r_moty-to-s_mask rotation matrix
function r2m_op(r_moty::Crystal, query::Crystal, m2r_isomorphism::Array{Int}, s_mask_atoms::Atoms)::Array{Float64,2}
    if m2r_isomorphism == Int[]
        return Matrix{Int}(I, 3, 3) # if no actual isom, skip OP and return identity
    end
    # r_moty subset in Cartesian
    A = r_moty.box.f_to_c * r_moty.atoms[m2r_isomorphism].coords.xf
    # s_mask in Cartesian
    B = query.box.f_to_c * s_mask_atoms.coords.xf
    # get rotation matrix
    return orthogonal_procrustes(A, B)
end


# Transforms r_moty according to two rotation matrices and a translational offset
function xform_r_moty(r_moty::Crystal, rot_r2m::Array{Float64,2}, rot_s2p::Array{Float64,2}, parent_offset::Array{Float64}, parent::Crystal)::Crystal
    # put r_moty into cartesian space
    atoms = Atoms{Cart}(length(r_moty.atoms.species), r_moty.atoms.species,
        Cart(r_moty.atoms.coords, r_moty.box))
    # transformation 1: rotate r_moty to align with query
    atoms.coords.x[:,:] = rot_r2m * atoms.coords.x
    # transformation 2: rotate to align with parent_subset
    atoms.coords.x[:,:] = rot_s2p * atoms.coords.x
    # transformation 3: translate to align with original parent center
    atoms.coords.x .+= parent.box.f_to_c * parent_offset
    # cast atoms back to Frac
    xrm = Crystal(r_moty.name, parent.box, Atoms{Frac}(length(atoms.species),
        atoms.species, Frac(atoms.coords, parent.box)), Charges{Frac}(0))
    # restore bonding network
    for e in edges(r_moty.bonds)
        add_edge!(xrm.bonds, src(e), dst(e))
    end
    return xrm
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
# of transformed r_moty's (xrms) along with the new fragments
function accumulate_bonds!(bonds::Array{Tuple{Int,Int}}, s2p_isom::Array{Int},
        parent::Crystal, m2r_isom::Union{Array{Int},Nothing}, xrm::Union{Crystal,Nothing}, count_xrms::Int)
    # skip bond accumulation for null replacement
    if m2r_isom == Int[] || isnothing(m2r_isom)
        return
    end
    # bonds between new fragments and parent
    # loop over s2p_isom
    for (s, p) in enumerate(s2p_isom)
        # in case the replacement moiety is smaller than the search moiety
        if s > length(m2r_isom)
            break
        end
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


# generates data for effecting a series of replacements
function build_replacement_data(configs::Array{Tuple{Int,Int}}, search::Search,
        parent::Crystal, query::Crystal, r_moty::Crystal, mask::Crystal, s′_in_r::Search)::Tuple{Array{Crystal},Array{Int},Array{Tuple{Int,Int}}}
    xrms = Crystal[]
    del_atoms = Int[]
    bonds = Tuple{Int,Int}[] # tuple (i,j) encodes a parent[i] -> xrms[k][j] bond
    # parent bonds
    for e in edges(parent.bonds) # loop over parent structure bonds
        push!(bonds, (src(e), dst(e))) # preserve them
    end
    # generate transformed replace moiety (xrm), ID which atoms to delete,
    # and accumulate list of bonds for each replacement configuration
    for config in configs
        # find isomorphism
        s2p_isom = search.results[config[1]].isomorphism[config[2]]
        # find parent subset
        parent_subset = deepcopy(parent[s2p_isom])
        # adjust coordinates for periodic boundaries
        adjust_for_pb!(parent_subset)
        # record the center of parent_subset so we can translate back later
        parent_subset_center = geometric_center(parent_subset)
        # shift to align centers at origin
        center_on_self!.([parent_subset, query])
        # orthog. Procrustes for query-to-parent and mask-to-replacement alignments
        rot_s2p = s2p_op(query, parent_subset)
        xrm = nothing
        m2r_isom = nothing
        if nb_isomorphisms(s′_in_r) ≠ 0
            # choose best r2m by evaluating RMSD for all possibilities
            rot_r2m_err = Inf
            for m2r_isom′ ∈ [s′_in_r.results[i].isomorphism[1] for i ∈ 1:nb_locations(s′_in_r)]
                # shift all r_moty nodes according to center of isomorphic subset
                r_moty′ = deepcopy(r_moty)
                r_moty′.atoms.coords.xf .-= geometric_center(r_moty[m2r_isom′])
                # get OP rotation matrix to align replacement onto query
                rot_r2m = r2m_op(r_moty, query, m2r_isom′, mask.atoms)
                # transform r_moty by rot_r2m, rot_s2p, and parent_subset_center, align
                # to parent (this is now a potential crystal to add)
                xrm′ = xform_r_moty(r_moty′, rot_r2m, rot_s2p, parent_subset_center, parent)
                # check error and keep xrm & m2r_isom if better than previous best error
                rot_r2m_err′ = rmsd(xrm′.atoms.coords.xf[:, m2r_isom′], mask.atoms.coords.xf[:, :])
                if rot_r2m_err′ < rot_r2m_err
                    m2r_isom = m2r_isom′
                    rot_r2m_err = rot_r2m_err′
                    xrm = xrm′
                end
            end
            # add optimal xrm to array
            push!(xrms, xrm)
        end
        # push obsolete atoms to array
        for x in s2p_isom
            push!(del_atoms, x) # this can create duplicates; remove them later
        end
        # clean up del_atoms
        del_atoms = unique(del_atoms)
        # accumulate bonds
        accumulate_bonds!(bonds, s2p_isom, parent, m2r_isom, xrm, length(xrms))
    end
    return xrms, del_atoms, bonds
end


## Search function (exposed)
@doc raw"""
    substructure_search(query, parent; disconnected_component=false)

Searches for a substructure within a `Crystal` and returns a `Search` struct
containing all identified subgraph isomorphisms.  Matches are made on the basis
of atomic species and chemical bonding networks, including bonds across unit cell
periodic boundaries.  The search moiety may optionally contain markup for
designating atoms to replace with other moieties.

# Arguments
- `query::Crystal` the search moiety
- `parent::Crystal` the parent structure
- `disconnected_component::Bool=false` if true, disables substructure searching and performs only exact matching
"""
function substructure_search(query::Crystal, parent::Crystal; disconnected_component::Bool=false)::Search
    # Make a copy w/o R tags for searching
    moty = deepcopy(query)
    untag_r_group!(moty)
    # Get array of configuration arrays
    configs = find_subgraph_isomorphisms(moty.bonds, moty.atoms.species, parent.bonds, parent.atoms.species, disconnected_component)
    df = DataFrame(p_subset=[sort(c) for c in configs], isomorphism=configs)
    locs = Int[]
    isoms = Array{Int}[]
    for (i, loc) in enumerate(groupby(df, :p_subset))
        for isom in loc.isomorphism
            push!(locs, i)
            push!(isoms, isom)
        end
    end
    results = groupby(DataFrame(location=locs, isomorphism=isoms), :location)
    return Search(SearchTerms(parent, query), results)
end


## Internal method for performing substructure replacements
function _substructure_replace(query::Crystal, r_moty::Crystal, parent::Crystal,
        search::Search, configs::Array{Tuple{Int,Int}},
        new_xtal_name::String)::Crystal
    # configs must all be unique
    @assert length(configs) == length(unique(configs)) "configs must be unique"
    # mutation guard
    query, r_moty = deepcopy.([query, r_moty])
    # if there are no replacements to be made, just return the parent
    if nb_isomorphisms(search) == 0
        @warn "No replacements to be made."
        return parent
    end
    # determine s_mask (which atoms from query are NOT in r_moty?)
    mask = query[idx_filter(query, r_group_indices(query))]
    # get isomrphism between query/mask and r_moty
    s′_in_r = mask ∈ r_moty
    if nb_isomorphisms(s′_in_r) ≠ 0
        m2r_isom = s′_in_r.results[1].isomorphism[1]
        # shift all r_moty nodes according to center of isomorphic subset
        r_moty.atoms.coords.xf .-= geometric_center(r_moty[m2r_isom])
    end
    # loop over configs to build replacement data
    xrms, del_atoms, bonds = build_replacement_data(configs, search, parent, query, r_moty, mask, s′_in_r)
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
    return new_xtal
end


## Find/replace function (exposed)
@doc raw"""
    substructure_replace(search, r_moty, nb_loc=2)

Inserts `r_moty` into a parent structure according to `search` and `kwargs`.
A valid replacement scheme must be selected by assigning one or more of the optional
`kwargs`.  Returns a new `Crystal` with the specified modifications (returns
`search.search.parent` if no replacements are made)

# Arguments
- `search::Search` the `Search` for a substructure moiety in the parent crystal
- `r_moty::Crystal` the moiety to use for replacement of the searched substructure
- `rand_all::Bool` set `true` to select random replacement at all matched locations
- `nb_loc::Int` assign a value to select random replacement at `nb_loc` random locations
- `loc::Array{Int}` assign value(s) to select specific locations for replacement.  If `ori` is not specified, replacement orientation is random.
- `ori::Array{Int}` assign value(s) when `loc` is assigned to specify exact configurations for replacement.
- `name::String` assign to give the generated `Crystal` a name ("new_xtal" by default)
"""
function substructure_replace(search::Search, r_moty::Crystal; rand_all::Bool=false,
        nb_loc::Int=0, loc::Array{Int}=Int[], ori::Array{Int}=Int[],
        name::String="new_xtal")::Crystal
    # handle input
    if rand_all # random replacement at each location
        nb_loc = nb_locations(search)
        loc = [1:nb_loc...]
        ori = [rand(1:nb_configs_at_loc(search)[i]) for i in loc]
    # random replacement at nb_loc random locations
    elseif nb_loc > 0 && ori == Int[] && loc == Int[]
        loc = sample([1:nb_locations(search)...], nb_loc, replace=false)
        ori = [rand(1:nb_configs_at_loc(search)[i]) for i in loc]
    elseif ori ≠ Int[] && loc ≠ Int[] # specific replacements
        @assert length(loc) == length(ori) "one orientation per location"
        nb_loc = length(ori)
    elseif loc ≠ Int[] # random replacement at specific locations
        nb_loc = length(loc)
        ori = [rand(1:nb_configs_at_loc(search)[i]) for i in loc]
    else
        @error "Invalid or missing replacement scheme."
    end
    # generate configuration tuples (location, orientation)
    configs = Tuple{Int,Int}[(loc[i], ori[i]) for i in 1:nb_loc]
    # process replacements
    return _substructure_replace(search.search.query, r_moty, search.search.parent,
        search, configs, name)
end

substructure_replace(search::Search, r_moty::Nothing; kwargs...) = substructure_replace(search, moiety(nothing); kwargs...)
