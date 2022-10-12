"""
    search = Search(parent, query, results)

Stores the `parent` and `query` used for a substructure search and the results (isomorphisms) of the subgraph matching algorithm.

## attributes

  - `search.parent::Crystal`                           # the parent in the search
  - `search.query::Crystal`                            # the query in the search
  - `search.isomorphisms::Vector{Vector{Vector{Int}}}` # the query-to-parent correspondences

The isomorphisms are grouped by location in the parent `Crystal` and can be examined using `nb_isomorphisms`, `nb_locations`, and `nb_ori_at_loc`.

Subgraph isomorphisms are encoded like

    isom = search.isomorphisms[i_loc][i_ori] = [7, 21, 9]

where `isom[k]` is the index of the atom in `search.parent` corresponding to atom `k` in `search.query` for the isomorphism at location `i_loc` and orientation `i_ori`.
"""
struct Search
    parent::Crystal
    query::Crystal
    isomorphisms::Vector{Vector{Dict{Int, Int}}}
end

function Base.show(io::IO, s::Search)
    begin
        println(io, s.query.name, " ∈ ", s.parent.name)
        print(io, nb_isomorphisms(s), " hits in ", nb_locations(s), " locations.")
    end
end

"""
    nb_isomorphisms(search::Search)

Returns the number of isomorphisms found in the specified `Search`

# Arguments

  - `search::Search` a substructure `Search` object
"""
function nb_isomorphisms(search::Search)::Int
    return sum(nb_ori_at_loc(search))
end

"""
    nb_locations(search::Search)

Returns the number of unique locations in the `parent` (sets of atoms in the `parent`) at which the
specified `Search` results contain isomorphisms.

# Arguments

  - `search::Search` a substructure `Search` object
"""
function nb_locations(search::Search)::Int
    return length(search.isomorphisms)
end

"""
    nb_ori_at_loc(search)

Returns a array containing the number of isomorphic configurations at a given
location (collection of atoms) for which the specified `Search` results
contain isomorphisms.

# Arguments

  - `search::Search` a substructure `Search` object
"""
function nb_ori_at_loc(search::Search)::Array{Int}
    return length.(search.isomorphisms)
end

# extension of infix `in` operator for syntactic sugar
# this allows all of the following:
#    s ∈ g                   →    find the moiety in the crystal
#    [s1, s2] .∈ [g]         →    find each moiety in a crystal
#    s .∈ [g1, g2]           →    find the moiety in each crystal
#    [s1, s2] .∈ [g1, g2]    →    find each moiety in each crystal
(∈)(s::Crystal, g::Crystal) = substructure_search(s, g)

"""
    iso_structs = isomorphic_substructures(s::Search)::Crystal

Returns a crystal consisting of the atoms of the `parent` involved in subgraph isomorphisms in the search `s`
"""
function isomorphic_substructures(s::Search)::Crystal
    return s.parent[reduce(
        vcat,
        collect.(values.([s.isomorphisms[i][1] for i in 1:nb_locations(s)]))
    )]
end

"""
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
function substructure_search(
    query::Crystal,
    parent::Crystal;
    disconnected_component::Bool=false
)::Search
    if parent.atoms.n == 0
        return
    end

    @assert ne(parent.bonds) > 0 "The parent structure must have bonds. Use `infer_bonds!(xtal, pbc)` to create them."
    # Make a copy w/o R tags for searching
    moty = deepcopy(query)
    untag_r_group!(moty)
    # Get array of configuration arrays
    configs = find_subgraph_isomorphisms(
        moty.bonds,
        moty.atoms.species,
        parent.bonds,
        parent.atoms.species,
        disconnected_component
    )
    df = DataFrame(; p_subset=[sort(c) for c in configs], isomorphism=configs)

    results = Vector{Dict{Int, Int}}[]
    for (i, df_loc) in enumerate(groupby(df, :p_subset))
        q2p_loc = Dict{Int, Int}[]
        for isomorphism in df_loc.isomorphism
            q2p = Dict([q => p for (q, p) in enumerate(isomorphism)])
            push!(q2p_loc, q2p)
        end
        push!(results, q2p_loc)
    end
    return Search(parent, query, results)
end
