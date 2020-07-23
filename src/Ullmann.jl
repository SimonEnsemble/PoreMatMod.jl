# ullmann.jl
# Adrian Henle, 2020

module Ullmann
export find_subgraph_isomorphisms
using PorousMaterials, LightGraphs

# "compatibility matrix"
# M₀[α, β] = 1 if any only if:
#     deg(β ∈ graph) ≥ deg(α ∈ subgraph)
#      and
#     species(α ∈ subgraph) == species(β ∈ graph)
function compatibility_matrix(subgraph::SimpleGraph, subgraph_species::Array{Symbol, 1},
                              graph::SimpleGraph,    graph_species::Array{Symbol, 1})::Array{Bool, 2}
    @debug "Finding M₀..."
    # allocate M. rows correspond to subgraph nodes, columns to graph nodes.
    M₀ = zeros(Bool, nv(subgraph), nv(graph))
    for α ∈ 1:nv(subgraph) # Loop over rows (subgraph nodes)
        for β ∈ 1:nv(graph) # Loop over columns (graph nodes)
            # Record Bool for each (i,j): true if atom species match and graph node degree is sufficient.
            if (degree(graph, β) ≥ degree(subgraph, α)) && (subgraph_species[α] == graph_species[β])
                M₀[α, β] = true # cannot rule out correspondence, so, true, for all we know at this point, these could correspond.
            end
        end
    end
    return M₀
end

# list of nodes β ∈ graph that could possibly correpond with node α ∈ subgraph
function candidate_list(M::Array{Bool, 2}, α::Int)::Array{Int, 1}
    return [β for β = 1:size(M, 2) if M[α, β]]
end

# does node α have possible candidate matches in the graph?
function has_candidates(M::Array{Bool, 2}, α::Int)::Bool
    for β ∈ 1:size(M, 2) # loop over graph nodes
        if M[α, β]
            return true
        end
    end
    # if made it this far, no graph node could possible correspond to α ∈ subgraph =/
    return false
end

function is_isomorphism(M::Array{Bool, 2})::Bool
    # (1) each row of M, corresponding to a node α ∈ subgraph, contains exactly one 1.
    #     i.e., every subgraph node has exactly one correspondence
    for α ∈ 1:size(M, 1)
        if sum(M[α, :]) != 1
            return false
        end
    end
    # (2) no column of M, corresponding to a node β ∈ graph, contains more than one 1.
    #     i.e., a graph node does not correspond to more than 1 subgraph node.
    for β ∈ 1:size(M, 2)
        if sum(M[:, β]) > 1
            return false
        end
    end
    # if made it this far, it is indeed an isomorphism.
    return true
end

# idea here:
#   if any subgraph node α has no possible correspondence w/ a node β in the graph, no point in continuing
#   return true iff M has no empty candidate lists for subgraph nodes.
function possibly_contains_isomorphism(M::Array{Bool, 2})::Bool
    @debug "Validating M: $(M)"
    for α ∈ 1:size(M, 1) # loop over subgraph nodes
        if ! has_candidates(M, α)
            return false # subgraph node α cannot be assigned! =-O no point in continuing.
        end
    end
    return true # M may be the intersection of one or more solutions.
end

function prune!(M::Array{Bool, 2}, subgraph::SimpleGraph, graph::SimpleGraph)
    pruned = true # to enter while loop
    while pruned
        pruned = false
        for α ∈ 1:size(M, 1) # loop thru subgraph nodes
            # get neighbors of node α
            neighbors_of_α = neighbors(subgraph, α)
            # loop thru candidate matches β ∈ graph for this subgraph node α
            for β ∈ candidate_list(M, α)
                neighbors_of_β = neighbors(graph, β)
                # now, suppose α ∈ subgraph and β ∈ graph correspond...
                for x ∈ neighbors_of_α
                    # if there is no neighbor of β that could correspond to x, neighbor of α, then, contradiction.
                    if length(intersect(candidate_list(M, x), neighbors_of_β)) == 0
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
    M[α, β] = true # assign correspondence at intersection
end

# soln:
#    soln[α ∈ subgraph] = β ∈ graph where α corresponds to β
function depth_first_search(α::Int, subgraph::SimpleGraph, graph::SimpleGraph, 
                            M::Array{Bool, 2}, soln::Array{Int, 1}, 
                            β_mapped::Array{Bool, 1}, solns::Array{Array{Int, 1}, 1})
    # if reached here from previous solution, exit.
    if α > size(M, 1)
        return nothing
    end

    # loop thru un-assigned graph nodes β that could possibly correspond to subnode α
    for β ∈ candidate_list(M, α)
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
            depth_first_search(α + 1, subgraph, graph, M′, soln, β_mapped, solns)
        end

        β_mapped[β] = false
        soln[α] = 0
    end
    return nothing
end

function find_subgraph_isomorphisms(subgraph::SimpleGraph, subgraph_species::Array{Symbol, 1},
                                    graph::SimpleGraph,    graph_species::Array{Symbol, 1})
    # store list of solutions here
    solns = Array{Array{Int, 1}, 1}()
    # encodes an isomorhism. maps α ∈ subgraph --> β ∈ graph
    soln = [0 for i = 1:nv(subgraph)]
    # tell us which β ∈ graph are mapped already.
    #   entry β true iff β mapped
    β_mapped = [false for i = 1:nv(graph)]
    # initial compatability matrix based on degrees of nodes and species
    M₀ = compatibility_matrix(subgraph, subgraph_species, graph, graph_species)
    depth_first_search(1, subgraph, graph, M₀, soln, β_mapped, solns)
    return solns
end

end # module
