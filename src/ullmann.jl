# PMUllmann.jl
# Adrian Henle, 2020

"""
    Ullmann algorithm code for PorousMaterials.jl
"""

## TODO
"""
    # fix ullmann_DFS()
    # sort M° on node degree in ullmann_bijections()
    # validate on sufficient unique candidates?
    # other accelerations/boosts?
    # ullmann tests
    # API → MOFun.jl
    # metadata dicts
"""
## #

module Ullmann


    using PorousMaterials, LightGraphs, DataFrames, LightGraphs.LinAlg


    """
        Generates the initial matrix based on node degrees.
        M will be a matrix indicating if the jth node of graph has
        sufficient degree to correspond with the ith node of subgraph.
    """
    function correspondence_matrix(subgraph::SimpleGraph, 𝒫s::DataFrame,
                                   graph::SimpleGraph, 𝒫g::DataFrame)::Array{Bool, 2}
        M = zeros(Bool, nv(subgraph), nv(graph))
        for i ∈ 1:nv(subgraph)
            for j ∈ 1:nv(graph)
                M[i, j] = 𝒫g.degree[j] ≥ 𝒫s.degree[i] && 𝒫g.species[j] == 𝒫s.species[i]
            end
        end
        return M
    end


    """
        Generates a new matrix from M0 by supposing x corresponds to y
    """
    function suppose_correspondence(M0::Array{Bool, 2},
                                    subgraph_node::Int,
                                    graph_node::Int)::Array{Bool, 2}
        M = M0 # deepcopy()?
        M[subgraph_node, :] .= false
        M[:, graph_node] .= false
        M[subgraph_node, graph_node] = true
        return M
    end


    """
        Gives the list of candidates for mapping the node according to M
        Argument `graph` must be a symbol, either :graph or :subgraph,
        to indicate to which graph the node belongs.
    """
    function candidate_list(M::Array{Bool,2},node_index::Int,graph::Symbol)::Array{Bool,2}
        # Wrong! Need to ONLY return the indexes of 1's, not the 1's and 0's
        if graph == :subgraph
            return M[node_index, :]
        elseif graph == :graph
            return M[:, node_index]
        else
            @error "Must provide either :graph or :subgraph."
        end
    end


    """
        Determines if M has no empty candidate lists for subgraph nodes.
    """
    function validate_correspondence_matrix(M::Array{Bool, 2})::Bool

## TODO add checking for sufficient unique destinations?

        for S_node ∈ 1:size(M, 1)
            if candidate_list(M, S_node, :subgraph) == []
                return false
            end
        end

        return true
    end


    """
        Determines if M is a 1-to-1 mapping solution
    """
    function is_solution(M::Array{Bool, 2})::Bool

        for i ∈ 1:size(M, 1)
            if sum(candidate_list(M, i, :subgraph)) ≠ 1
                return false
            end
        end

        for i ∈ 1:size(M, 2)
            if sum(candidate_list(M, i, :graph)) > 1
                return false
            end
        end

        return true
    end


    """
        Returns list of neighbors of input node from its graph's adjacency matrix.
    """
    function neighbors(w::Int, A::Array{Bool, 2})::Array{Int}
        return [v for v ∈ 1:size(A, 1) if A[v, w]]
    end


    """
        Performs Ullman algorithm refinement until M is stable
    """
    function refine_correspondence_matrix(  As::Array{Bool, 2},
                                            Ag::Array{Bool, 2},
                                            M0::Array{Bool, 2} )::Array{Bool, 2}

        M = M0 # don't need copying here! Make this a modify!() function

        # Wrong!  The core proposition is not implemented right.
        for y ∈ size(M0, 1)
            for x ∈ candidate_list(M, y, :graph)
                for z ∈ neighbors(y, As)
                    if candidate_list(M, z, :subgraph)
                        M[y, x] = 0
                    end
                end
            end
        end

        if M ≠ M0 # Expensive!! Change to Bool flag at line 159
            M = refine_correspondence_matrix(As, Ag, M)
        end

        return M
    end


    """
        Performs depth-first search for Ullmann's algorithm.
    """
    function ullmann_DFS(M0::Array{Bool, 2}, AS::Array{Bool, 2},
                         AG::Array{Bool, 2})::Array{Array{Bool, 2}}
        ℳ = []
        x1 = -1
        y1 = -1

        # Find first possible but not necessary correspondence
## TODO
# Wrong. This always finds upper-leftmost "hit" instead of next correct one.
        for y ∈ 1:size(M0, 1)
            for x ∈ 1:size(M0, 2)
                if M0[y, x]
                    x1 = x
                    break
                end
                if x1 ≠ -1
                    y1 = y
                    break
                end
            end
        end

        if !(validate_correspondence_matrix(M0))
            return []
        else
            M′ = suppose_correspondence(M0, y1, x1)
            M′ = refine_correspondence_matrix(As, Ag, M′)

            if is_solution(M′)
                return append!(ℳ, M′)
            else
                return append!(ℳ, ullmann_DFS(M′, AS, AG))
            end
        end
    end


    """
        Returns all correspondence matrix bijections of subgraph to a subset of graph.
    """
    function ullmann_bijections(subgraph::SimpleGraph,
                                subgraph_species::Array{Symbol},
                                graph::SimpleGraph,
                                graph_species::Array{Symbol})::Array{Array{Bool,2}}
        # Build metadata dictionaries
        𝒫s = DataFrame(index = 1:nv(subgraph), species = subgraph_species, degree = degree(subgraph))
        𝒫g = DataFrame(index = 1:nv(graph), species = graph_species, degree = degree(graph))
        # Get initial candidate correspondence matrix and adjacency matrices
        M° = correspondence_matrix(subgraph, 𝒫s, graph, 𝒫g)
        𝒜s = Array{Bool}(LinAlg.adjacency_matrix(subgraph))
        𝒜g = Array{Bool}(LinAlg.adjacency_matrix(graph))
## TODO node degree sorting descending on S
        # Perform depth-first search
        ℳ = ullmann_DFS(M°, 𝒜s, 𝒜g)
        # Return solutions
        return ℳ
    end

    ## PMUllmann API




    """
    	Returns a Crystal wherein all instances of search_moiety from parent_structure
    	are replaced by replace_moiety.  Uses Ullmann's algorithm for matching.
    """
    function subgraph_find_replace(search_moiety::Crystal,
    							   replace_moiety::Crystal,
    							   parent_structure::Crystal)::Array{Array{Bool, 2}}
    ## TODO implement
    	matches = subgraph_isomorphisms(search_moiety, parent_structure)
    	return nothing
    end
    """
    	Returns all bijections that map the bonding network from search_moiety onto
    	parent_structure via Ullmann's algorithm
    """
    function subgraph_isomorphisms(search_moiety::Crystal,
    							   parent_structure::Crystal)::Array{Array{Bool, 2}}
    	return ullmann_bijections(search_moiety.bonds, search_moiety.atoms.species,
    							  parent_structure.bonds, parent_structure.atoms.species)
    end

export subgraph_isomorphisms

end
