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


# extension of infix `in` operator for expressive searching
# this allows all of the following:
#    s ∈ g                   →    find the moiety in the crystal
#    [s1, s2] .∈ [g]         →    find each moiety in a crystal
#    s .∈ [g1, g2]           →    find the moiety in each crystal
#    [s1, s2] .∈ [g1, g2]    →    find each moiety in each crystal
(∈)(s::Crystal, g::Crystal) = substructure_search(s, g)


"""
`iso_structs = isomorphic_substructures(s::Search)::Crystal`

Returns a crystal consisting of the atoms involved in subgraph isomorphisms in the search `s`
"""
function isomorphic_substructures(s::Search)::Crystal
    return s.search.parent[reduce(vcat, [s.results[i].isomorphism[1] for i in 1:nb_locations(s)])]
end


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
