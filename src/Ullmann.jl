# "compatibility matrix"
# M₀[α, β] = 1 if any only if:
#     deg(β ∈ graph) ≥ deg(α ∈ subgraph)
#      and
#     species(α ∈ subgraph) == species(β ∈ graph)
function compatibility_matrix(
    subgraph::MetaGraph,
    subgraph_species::Array{Symbol, 1},
    graph::MetaGraph,
    graph_species::Array{Symbol, 1},
    disconnected_component::Bool
)::Array{Bool, 2}
    # allocate M. rows correspond to subgraph nodes, columns to graph nodes.
    M₀ = zeros(Bool, nv(subgraph), nv(graph))

    # Get adjacency matrices and 4th/8th degree path matrices
    adjmat_S = adjacency_matrix(subgraph)
    adjmat_G = adjacency_matrix(graph)
    deg_S = adjmat_S^2
    deg_G = adjmat_G^2
    path4_S = deg_S^2
    path4_G = deg_G^2
    path8_S = path4_S^2
    path8_G = path4_G^2

    if !disconnected_component # search for substructures
        @inbounds for β in 1:nv(graph) # Loop over rows (subgraph nodes)
            @inbounds for α in 1:nv(subgraph) # Loop over columns (graph nodes)
                # Record Bool for each (i,j): true if atom species match, graph node degree is sufficient, and 4 and 8 length self-paths are sufficient.
                M₀[α, β] =
                    subgraph_species[α] == graph_species[β] &&
                    deg_G[β, β] ≥ deg_S[α, α] &&
                    path4_G[β, β] ≥ path4_S[α, α] &&
                    path8_G[β, β] ≥ path8_S[α, α]
            end
        end
    else # search only for exact, isolated matches (no substructures)
        @inbounds for β in 1:nv(graph) # Loop over rows (subgraph nodes)
            @inbounds for α in 1:nv(subgraph) # Loop over columns (graph nodes)
                # Record Bool for each (i,j): true if atom species match, and graph node degree matches.
                M₀[α, β] =
                    subgraph_species[α] == graph_species[β] && deg_G[β, β] == deg_S[α, α]
            end
        end
    end
    return M₀
end

# list of nodes β ∈ graph that could possibly correpond with node α ∈ subgraph
function candidate_list(M::Array{Bool, 2}, α::Int)::Array{Int, 1}
    @inbounds @views return findall(M[α, :])
end

# does node α have possible candidate matches in the graph?
function has_candidates(M::Array{Bool, 2}, α::Int)::Bool
    @inbounds return any([M[α, β] for β in 1:size(M, 2)])
end

function is_isomorphism(M::Array{Bool, 2})::Bool
    # (1) each row of M, corresponding to a node α ∈ subgraph, contains exactly one 1.
    #     i.e., every subgraph node has exactly one correspondence
    # (2) no column of M, corresponding to a node β ∈ graph, contains more than one 1.
    #     i.e., a graph node does not correspond to more than 1 subgraph node.
    @inbounds return !(
        any([sum(M[α, :]) ≠ 1 for α in 1:size(M, 1)]) || any([sum(M[:, β]) > 1 for β in 1:size(M, 2)])
    )
end

# idea here:
#   if any subgraph node α has no possible correspondence w/ a node β in the graph, no point in continuing
#   return true iff M has no empty candidate lists for subgraph nodes.
function possibly_contains_isomorphism(M::Array{Bool, 2})::Bool
    @inbounds return all([has_candidates(M, α) for α in 1:size(M, 1)])
end

function prune!(M::Array{Bool, 2}, subgraph::MetaGraph, graph::MetaGraph)
    pruned = true # to enter while loop
    while pruned
        pruned = false
        @inbounds for α in 1:size(M, 1) # loop thru subgraph nodes
            # get neighbors of node α
            neighbors_of_α = neighbors(subgraph, α)
            # loop thru candidate matches β ∈ graph for this subgraph node α
            @inbounds for β in candidate_list(M, α)
                neighbors_of_β = neighbors(graph, β)
                # now, suppose α ∈ subgraph and β ∈ graph correspond...
                @inbounds for x in neighbors_of_α
                    # if there is no neighbor of β that could correspond to x, neighbor of α, then, contradiction.
                    if isdisjoint(candidate_list(M, x), neighbors_of_β)
                        M[α, β] = false
                        pruned = true
                    end
                end
            end
        end
    end
end

function assign_correspondence!(M::Array{Bool, 2}, α::Int, β::Int)
    M[α, :] .= false # zero out row of subgraph node
    M[:, β] .= false # zero out column of graph node
    return M[α, β] = true # assign correspondence at intersection
end

# soln:
#    soln[α ∈ subgraph] = β ∈ graph where α corresponds to β
function depth_first_search!(
    α::Int,
    subgraph::MetaGraph,
    graph::MetaGraph,
    M::Array{Bool, 2},
    soln::Array{Int, 1},
    β_mapped::Array{Bool, 1},
    solns::Array{Array{Int, 1}, 1}
)
    # if reached here from previous solution, exit.
    if α > size(M, 1)
        return nothing
    end

    # loop thru un-assigned graph nodes β that could possibly correspond to subnode α
    @inbounds for β in candidate_list(M, α)
        # if βraph is already mapped, not a viable solution.
        #   (not sure if necessary, i.e. if M already knows this)
        if β_mapped[β]
            continue
        end

        # make a copy so we can restore later.
        M′ = deepcopy(M)

        # explore scenario where α ∈ subgraph corresponds to β ∈ graph
        assign_correspondence!(M′, α, β)
        soln[α] = β
        β_mapped[β] = true

        # prune tree
        prune!(M′, subgraph, graph)

        # if we reached bottom of tree, iso-morphism is found!
        if α == size(M′, 1)
            # do we hv to check if it's a sol'n or is it guarenteed? why prune then?
            if is_isomorphism(M′)
                push!(solns, deepcopy(soln))
            end
            # don't return b/c we need to look at other candidates
        end

        if M′[α, β] && possibly_contains_isomorphism(M′)
            # we've assigned α, go deeper in the depth first search
            depth_first_search!(α + 1, subgraph, graph, M′, soln, β_mapped, solns)
        end

        β_mapped[β] = false
        soln[α] = 0
    end
end

@doc raw"""
returns an array of arrays, each containing one unique subgraph isomorphism
"""
function find_subgraph_isomorphisms(
    subgraph::MetaGraph,
    subgraph_species::Array{Symbol, 1},
    graph::MetaGraph,
    graph_species::Array{Symbol, 1},
    disconnected_component::Bool=false
)
    # store list of solutions here
    solns = Array{Array{Int, 1}, 1}()
    # encodes an isomorhism. maps α ∈ subgraph --> β ∈ graph
    soln = [0 for _ in 1:nv(subgraph)]
    # tell us which β ∈ graph are mapped already.
    #   entry β true iff β mapped
    β_mapped = [false for _ in 1:nv(graph)]
    # initial compatability matrix based on degrees of nodes and species
    M₀ = compatibility_matrix(
        subgraph,
        subgraph_species,
        graph,
        graph_species,
        disconnected_component
    )
    depth_first_search!(1, subgraph, graph, M₀, soln, β_mapped, solns)
    return solns
end
