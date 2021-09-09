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
`iso_structs = isomorphic_substructures(s::Search)::Crystal`

Returns a crystal consisting of the atoms involved in subgraph isomorphisms in the search `s`
"""
function isomorphic_substructures(s::Search)::Crystal
    return s.search.parent[reduce(vcat, [s.results[i].isomorphism[1] for i in 1:nb_locations(s)])]
end


# Translates all atoms in xtal such that xtal[1] is in its original position
# and the rest of xtal is in its nearest-image position relative to xtal[1]
function adjust_for_pb!(xf::Matrix{Float64}; xtal_name::String="")
    # record position vector of xtal[1]
    origin_offset = deepcopy(xf[:, 1])
    # loop over atom indices and correct coordinates
    for i in 1:size(xf, 2)
        # move atoms near the origin for nearest-image calculation
        dxf = xf[:, i] .- origin_offset
        # nearest_image! expects points to be within same or adjacent unit cells
        @assert all(abs.(dxf) .< 2) "Invalid xf coords $xtal_name"
        # resolve periodic boundaries (other vectors unchanged)
        nearest_image!(dxf)
        # return atoms to their [nearest-image] original positions
        xf[:, i] = dxf .+ origin_offset
    end
end

adjust_for_pb!(xtal::Crystal) = adjust_for_pb!(xtal.atoms.coords.xf, xtal_name=xtal.name)


# Gets the rotation matrix for aligning the replacement moiety onto a subset (isomorphic to masked query) of parent atoms
function r2p_op(replacement::Crystal, parent::Crystal, r2p_isom::Dict{Int,Int}, parent_subset_center::Vector{Float64})
    @assert replacement.atoms.n ≥ 3 && parent.atoms.n ≥ 3 "Parent and replacement must each be at least 3 atoms for SVD alignment."
    # get matrix A (replacement fragment coordinates)
    A = replacement.box.f_to_c * replacement.atoms[[r for (r,p) in r2p_isom]].coords.xf
    # prepare parent subset
    parent_subset = deepcopy(parent.atoms[[p for (r,p) in r2p_isom]].coords.xf)
    adjust_for_pb!(parent_subset)
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


# returns an Array containing the indices
function idx_filter(xtal::Crystal, subset::Array{Int})::Array{Int,1}
    return [i for i in 1:xtal.atoms.n if !(i ∈ subset)]
end


# tracks which bonds need to be made between the parent and array
# of transformed replacement's (xrms) along with the new fragments
function accumulate_bonds!(bonds::Array{Tuple{Int,Int}}, q2p_isom::Array{Int},
        parent::Crystal, m2r_isom::Union{Array{Int},Nothing}, xrm::Union{Crystal,Nothing}, count_xrms::Int)
    # skip bond accumulation for null replacement
    if m2r_isom == Int[] || isnothing(m2r_isom)
        return
    end
    # bonds between new fragments and parent
    # loop over q2p_isom
    for (s, p) in enumerate(q2p_isom)
        # in case the replacement moiety is smaller than the search moiety
        if s > length(m2r_isom)
            break
        end
        # find neighbors of parent_subset atoms
        n = LightGraphs.neighbors(parent.bonds, p)
        # loop over neighbors
        for nᵢ in n
            # if neighbor not in q2p_isom, must bond it to replacement replacement
            # of parent_subset atom
            if !(nᵢ ∈ q2p_isom)
                # ID the atom in replacement
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


# chains (m2r)⁻¹ -> m2q -> q2p to obtain r2p mapping
function map_r′_to_p(m2r, q2p, m2q)
    # get inverse of m2r_isom, i.e. r2m
    r2m = Dict([r => m for (m,r) in enumerate(m2r)])
    # chain r2m -> m2q
    r2q = Dict([r => m2q[m] for (r,m) in r2m])
    # chain r2q -> q2p
    r2p = Dict([r => q2p[q] for (r,q) in r2q])
    return r2p
end


# generates data for effecting a series of replacements
function build_replacement_data(configs::Vector{Tuple{Int,Int}}, q_in_p::Search,
        parent::Crystal, replacement::Crystal, q′_in_r::Search, m2q_key::Vector{Int}
        )::Tuple{Vector{Crystal},Vector{Int},Vector{Tuple{Int,Int}}}
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
        ori = config[2]
        if ori == 0
            ori = [1:nb_configs_at_loc(q_in_p)[loc]...]
        else
            ori = [ori]
        end
        xrm = nothing
        m2r_isom = nothing
        alignment_err = Inf
        q2p_isom = nothing
        for ori′ in ori
            # pull up specific isomorphism from query to parent
            q2p_isom′ = q_in_p.results[loc].isomorphism[ori′]
            # determine isomorphism from masked query to parent
            m2p_isom = q2p_isom′[m2q_key]
            # find parent subset
            parent_subset = deepcopy(parent[m2p_isom])
            # adjust coordinates for periodic boundaries
            adjust_for_pb!(parent_subset)
            # record the center of parent_subset so we can translate back later
            parent_subset_center = geometric_center(parent_subset)
            
            # orthog. Procrustes
            if nb_isomorphisms(q′_in_r) ≠ 0
                # choose best r2p by evaluating RMSD for all possibilities
                for i ∈ 1:nb_locations(q′_in_r)
                    for j ∈ 1:nb_configs_at_loc(q′_in_r)[i]
                        m2r_isom′ = q′_in_r.results[i].isomorphism[j]
                        # shift all replacement nodes according to center of isomorphic subset
                        replacement′ = deepcopy(replacement)
                        replacement′.atoms.coords.xf .-= geometric_center(replacement[m2r_isom′])
                        # determine mapping r2p
                        r2p_isom = map_r′_to_p(m2r_isom′, q2p_isom′, m2q_key)
                        # get OP rotation matrix to align replacement onto parent
                        rot_r2p = r2p_op(replacement, parent, r2p_isom, parent_subset_center)
                        # transform replacement by rot_r2p, and parent_subset_center (this is now a potential crystal to add)
                        xrm′ = xform_replacement(replacement′, rot_r2p, parent_subset_center, parent)
                        # check error and keep xrm & m2r_isom if better than previous best error
                        alignment_err′ = rmsd(xrm′.atoms.coords.xf[:, m2r_isom′], parent.atoms.coords.xf[:, m2q_key])
                        if alignment_err′ < alignment_err
                            m2r_isom = m2r_isom′
                            alignment_err = alignment_err′
                            xrm = xrm′
                            q2p_isom = q2p_isom′
                        end
                    end
                end
            else
                q2p_isom = q2p_isom′
            end
        end
        # add optimal xrm to array
        if !isnothing(xrm)
            push!(xrms, xrm)
        end
        # push obsolete atoms to array
        for x in q2p_isom
            push!(del_atoms, x) # this can create duplicates; remove them later
        end
        # clean up del_atoms
        del_atoms = unique(del_atoms)
        # accumulate bonds
        accumulate_bonds!(bonds, q2p_isom, parent, m2r_isom, xrm, length(xrms))
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
- `disconnected_component::Bool=false` if true, disables substructure searching (e.g. for finding guest molecules)
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
function _substructure_replace(q_in_p::Search, replacement::Crystal, configs::Array{Tuple{Int,Int}}, new_xtal_name::String)::Crystal
    query = q_in_p.search.query
    parent = q_in_p.search.parent
    # configs must all be unique
    @assert length(configs) == length(unique(configs)) "configs must be unique"
    # mutation guard
    query, replacement = deepcopy.([query, replacement])
    # if there are no replacements to be made, just return the parent
    if nb_isomorphisms(q_in_p) == 0
        @warn "No replacements to be made."
        return parent
    end
    # which atoms from query are in replacement?
    m2q_key = idx_filter(query, r_group_indices(query))
    # get isomrphism between query/mask and replacement
    q′_in_r = query[m2q_key] ∈ replacement
    # loop over configs to build replacement data
    xrms, del_atoms, bonds = build_replacement_data(configs, q_in_p, parent, replacement, q′_in_r, m2q_key)
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
    substructure_replace(search, replacement, nb_loc=2)

Inserts `replacement` into a parent structure according to `search` and `kwargs`.
Provide a replacement style `kwarg` to direct the location and orientation of replacements.
Default behavior is to seek the replacement operation with lowest RMSD on spatial alignment.
Returns a new `Crystal` with the specified modifications (returns `search.search.parent` if no replacements are made)

# Arguments
- `search::Search` the `Search` for a substructure moiety in the parent crystal
- `replacement::Crystal` the moiety to use for replacement of the searched substructure
- `random::Bool` set `true` to select random replacement orientations
- `nb_loc::Int` assign a value to select random replacement at `nb_loc` random locations
- `loc::Array{Int}` assign value(s) to select specific locations for replacement.  If `ori` is not specified, replacement orientation is random.
- `ori::Array{Int}` assign value(s) when `loc` is assigned to specify exact configurations for replacement. `0` values mean the configuration at that location should be selected for optimal alignment with the parent.
- `name::String` assign to give the generated `Crystal` a name ("new_xtal" by default)
- `verbose::Bool` set `true` to print console messages about the replacement(s) being performed
"""
function substructure_replace(search::Search, replacement::Crystal; random::Bool=false,
        nb_loc::Int=0, loc::Array{Int}=Int[], ori::Array{Int}=Int[], name::String="new_xtal", verbose::Bool=false)::Crystal
    # replacement at all locations (default)
    if nb_loc == 0 && loc == Int[] && ori == Int[]
        nb_loc = nb_locations(search)
        loc = [1:nb_loc...]
        if random
            ori = [rand(1:nb_configs_at_loc(search)[i]) for i in loc]
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
            ori = [rand(1:nb_configs_at_loc(search)[i]) for i in loc]
            if verbose
                @info "Replacing" q_in_p=search r=replacement.name mode="random ori @ $nb_loc loc"
            end
        else
            ori = zeros(Int, nb_loc)
            if verbose
                @info "Replacing" q_in_p=search r=replacement.name mode="optimal ori @ $nb_loc loc"
            end
        end
    # replacement at specific locations
    elseif loc ≠ Int[]
        nb_loc = length(loc)
        if random
            ori = [rand(1:nb_configs_at_loc(search)[i]) for i in loc]
            if verbose
                @info "Replacing" q_in_p=search r=replacement.name mode="random ori @ loc: $loc"
            end
        else
            ori = zeros(Int, nb_loc)
            if verbose
                @info "Replacing" q_in_p=search r=replacement.name mode="optimal ori @ loc: $loc"
            end
        end
    # specific replacements
    elseif ori ≠ Int[] && loc ≠ Int[]
        @assert length(loc) == length(ori) "one orientation per location"
        nb_loc = length(ori)
        if verbose
            @info "Replacing" q_in_p=search r=replacement.name mode="loc: $loc\tori: $ori"
        end
    end

    # generate configuration tuples (location, orientation)
    configs = Tuple{Int,Int}[(loc[i], ori[i]) for i in 1:nb_loc]
    # process replacements
    return _substructure_replace(search, replacement, configs, name)
end

substructure_replace(search::Search, replacement::Nothing; kwargs...) = substructure_replace(search, moiety(nothing); kwargs...)
