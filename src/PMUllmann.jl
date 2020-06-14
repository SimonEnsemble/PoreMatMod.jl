# PMUllmann.jl
# Adrian Henle, 2020

"""
    Ullmann algorithm code for PorousMaterials.jl
"""

## TODO
"""
    # fix ullmann_bijections()
    # fix ullmann_DFS()
    # sort M° on node degree in ullmann_bijections()
    # validate on sufficient unique candidates?
    # ullmann tests
    # API → MOFun.jl
    # add path_to_moieties to PorousMaterials
    # moiety() tests:
        # generate vtk
        # generate manual graph representation
        # verify that the two are identical
    # metagraphs
"""
## #

module PMUllmann


    using LightGraphs


    export ullmann_bijections


    """
        Make the module flexible for updating to MetaGraphs
    """
## TODO metagraphs.
    abstract type Graph <: SimpleGraph end


    """
        Returns a tuple of diagonal matrices (𝒜s, 𝒜g), the adjacency matrices
        for input subgraph and graph, respectively.
    """
    function adjacency_matrices(subgraph::Graph, graph::Graph)::Array{Bool, 2}
        return (LinAlg.adjacency_matrix(subgraph), LinAlg.adjacency_matrix(graph))
    end


    """
        Generates the initial matrix based on node degrees.
    """
    function correspondence_matrix(subgraph::Graph, graph::Graph)::Array{Bool, 2}

        M = zeros(nv(subgraph), nv(graph))
        Ds = degree(subgraph)
        Dg = degree(graph)

        for i ∈ 1:nv(subgraph)
            for j ∈ 1:nv(graph)
                M[i, j] = (Dg[j] - Ds[i] ≥ 0)
            end
        end

        return M
    end


    """
        Generates a new matrix from M0 by supposing x corresponds to y
    """
    function suppose_correspondence(    M0::Array{Bool, 2}},
                                        subgraph_node::Int,
                                        graph_node::Int )::Array{Bool, 2}

        M = M0
        M[subgraph_node, :] = 0
        M[:, graph_node] = 0
        M[subgraph_node, graph_node] = 1

        return M
    end


    """
        Gives the list of candidates for mapping the node according to M
        Argument `graph` must be a symbol, either :graph or :subgraph,
        to indicate to which graph the node belongs.
    """
    function candidate_list(    M::Array{Bool, 2},
                                node_index::Int,
                                graph::Symbol )::Array{Bool, 2}

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
    function is_solution(M::Array{Bool, 2}})::Bool

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

        N = []
        for v ∈ 1:size(A, 1)
            if A[v] == 1
                append!(N, A[v])
            end
        end

        return N
    end


    """
        Performs Ullman algorithm refinement until M is stable
    """
    function refine_correspondence_matrix(  As::Array{Bool, 2},
                                            Ag::Array{Bool, 2},
                                            M0::Array{Bool, 2} )::Array{Bool, 2}

        M = M0

        for y ∈ size(M0, 1)
            for x ∈ candidate_list(M, y, :graph)
                for z ∈ neighbors(y, As)
                    if candidate_list(M, z, :subgraph)
                        M[y, x] = 0
                    end
                end
            end
        end

        if M ≠ M0
            M = refine_correspondence_matrix(As, Ag, M)
        end

        return M
    end


    """
        Performs depth-first search for Ullmann's algorithm.
    """
    function ullmann_DFS(   M0::Array{Bool, 2},
                            AS::Array{Bool, 2},
                            AG::Array{Bool, 2} )::Array{Array{Bool, 2}}

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

        if ¬(validate_correspondence_matrix(M0))
            return []
        else
            M′ = suppose_correspondence(M0, x1, y1)
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
    function ullmann_bijections(subgraph::Graph, graph::Graph)
## TODO refactor to take Crystal objects (will need to extract graphs)

        # Get initial candidate correspondence matrix and adjacency matrices
        M° = correspondence_matrix(subgraph, graph)
        (𝒜s, 𝒜g) = adjacency_matrices(subgraph, graph)
## TODO node degree sorting

        # Return set of subgraph isomorphic bijections via depth-first search
        return ullmann_DFS(M°, 𝒜s, 𝒜g)
    end
end
