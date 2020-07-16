# ullmann.jl
# Adrian Henle, 2020

@doc raw"""
    Ullmann algorithm code for PorousMaterials.jl
"""
module Ullmann

    using PorousMaterials, LightGraphs, DataFrames, LightGraphs.LinAlg


    @doc raw"""
    function correspondence_matrix(subgraph::SimpleGraph, 𝒫s::DataFrame,
                                   graph::SimpleGraph, 𝒫g::DataFrame)::Array{Bool, 2}
        Generates the initial matrix based on node degrees.
        M will be a matrix indicating if the jth node of graph has
        sufficient degree to correspond with the ith node of subgraph.
    """
    function correspondence_matrix(subgraph::SimpleGraph, 𝒫s::DataFrame,
                                   graph::SimpleGraph, 𝒫g::DataFrame)::Array{Bool, 2}
        @debug "Finding M..."
        M = zeros(Bool, nv(subgraph), nv(graph))
        for i ∈ 1:nv(subgraph)
            for j ∈ 1:nv(graph)
@debug "($i,$j) $(𝒫g.degree[j]) ≥ $(𝒫s.degree[i]) && $(𝒫g.species[j]) == $(𝒫s.species[i]) : $(𝒫g.degree[j] ≥ 𝒫s.degree[i] && 𝒫g.species[j] == 𝒫s.species[i]))"
                M[i, j] = 𝒫g.degree[j] ≥ 𝒫s.degree[i] && 𝒫g.species[j] == 𝒫s.species[i]
            end
        end
        @debug "M: $(M)"
        return M
    end


    @doc raw"""
        function suppose_correspondence(M0::Array{Bool, 2},
                                    subgraph_node::Int,
                                    graph_node::Int)::Array{Bool, 2}
        Generates a new matrix from M0 by supposing x corresponds to y
    """
    function suppose_correspondence(M0::Array{Bool, 2},
                                    subgraph_node::Int,
                                    graph_node::Int)::Array{Bool, 2}
        @debug "Making supposition A[$(subgraph_node)]↔B[$(graph_node)]"
        M = deepcopy(M0)
        M[subgraph_node, :] .= false
        M[:, graph_node] .= false
        M[subgraph_node, graph_node] = true

        return M
    end


    @doc raw"""
        function candidate_list(M::Array{Bool,2},node_index::Int,graph::Symbol)::Array{Int}
        Gives the list of candidates for mapping the node according to M
        Argument `graph` must be a symbol, either :graph or :subgraph,
        to indicate to which graph the node belongs.
    """
    function candidate_list(M::Array{Bool,2},node_index::Int,graph::Symbol)::Array{Int}
        @debug "Getting candidate list for node $(node_index) of $(graph)"
        if graph == :subgraph
            logicals = M[node_index, :]
        elseif graph == :graph
            logicals = M[:, node_index]
        else
            @error "Must provide either :graph or :subgraph."
        end
        @debug "Logicals: $(logicals)"
        return [index for (index, logical) ∈ enumerate(logicals) if logical]
    end


    @doc raw"""
        function validate_M(M::Array{Bool, 2})::Bool
        Determines if M has no empty candidate lists for subgraph nodes.
    """
    function validate_M(M::Array{Bool, 2})::Bool
        @debug "Validating M: $(M)"
        for S_node ∈ 1:size(M, 1)
            if length(candidate_list(M, S_node, :subgraph)) == 0
                return false
            end
        end
        return true
    end


    @doc raw"""
        function is_solution(M::Array{Bool, 2})::Bool
        Determines if M is a 1-to-1 mapping solution
    """
    function is_solution(M::Array{Bool, 2})::Bool
        @debug "Testing for solution: $(M)"
        for i ∈ 1:size(M, 1)
            if length(candidate_list(M, i, :subgraph)) ≠ 1
                @debug "Not a solution."
                return false
            end
        end
        for i ∈ 1:size(M, 2)
            if length(candidate_list(M, i, :graph)) > 1 # reachable?
                @debug "Not a solution."
                return false
            end
        end
        @debug "Found a solution: $(M)"
        return true
    end


    @doc raw"""
        function neighbors(w::Int, A::Array{Bool, 2})::Array{Int}
        Returns list of neighbors of input node from its graph's adjacency matrix.
    """
    function neighbors(w::Int, A::Array{Bool, 2})::Array{Int}
        @debug "Finding neighbors of node $(w)."
        return [v for v ∈ 1:size(A, 1) if A[v, w]]
    end


    @doc raw"""
        function refine_M!(As::Array{Bool, 2}, Ag::Array{Bool, 2}, M::Array{Bool, 2})
        Performs Ullman algorithm refinement until M is stable
    """
    function refine_M!(As::Array{Bool, 2}, Ag::Array{Bool, 2}, M::Array{Bool, 2})
        @debug "Refining M: $(M)"
        M_altered = false
        for y ∈ size(M, 1) # Nodes in subgraph
            for x ∈ candidate_list(M, y, :subgraph)
                for z ∈ neighbors(y, As)
                    if length(intersect(candidate_list(M, z, :subgraph), neighbors(x, Ag))) == 0
                        M[y, x] = 0
                        @debug "M altered at [$(y), $(x)]."
                        M_altered = true
                    end
                end
            end
        end
        M = M_altered ? refine_M!(As, Ag, M) : M
    end


    @doc raw"""
        function ullmann_DFS(M::Array{Bool, 2}, AS::Array{Bool, 2}, AG::Array{Bool, 2},
                       i°::Int=1, j°::Int=1)::Union{Array{Array{Bool, 2}}, Nothing}
        Performs depth-first search for Ullmann's algorithm.
    """
    function ullmann_DFS(M::Array{Bool, 2}, AS::Array{Bool, 2}, AG::Array{Bool, 2},
                       i°::Int=1, j°::Int=1)::Dict{Int, Array{Bool, 2}}
        ℳ = Dict{Int, Array{Bool, 2}}()
        if validate_M(M)
            for i ∈ i°:size(M, 1) # Looping over subgraph nodes
                @debug "i: $i"
                if length(candidate_list(M, i, :graph)) > 1 # Only deal with columns on this row if it is not already solved
                    for j ∈ j°:size(M, 2) # Looping over graph nodes
                        @debug "j: $j"
                        if M[i, j]
                            M′ = suppose_correspondence(M, i, j)
                            refine_M!(AS, AG, M′)
                            if is_solution(M′)
                                @debug "Appending solution to ℳ: $(M′)"
                                ℳ[length(ℳ)+1] = copy(M′)
                                @debug "Solutions:" ℳ
                            elseif validate_M(M′)
                                @debug "Continuing search..."
                                ℳ[length(ℳ)+1] = ullmann_DFS(M′, AS, AG, i, j)
                            else
                                @debug "Reached leaf node."
                            end
                        end
                    end
                end
            end
        end
        return ℳ
    end


    @doc raw"""
        function ullmann_bijections(subgraph::SimpleGraph,
                                subgraph_species::Array{Symbol},
                                graph::SimpleGraph,
                                graph_species::Array{Symbol})::Union{Array{Array{Bool,2}}, Nothing}
        Returns all correspondence matrix bijections of subgraph to a subset of graph.
    """
    function ullmann_bijections(subgraph::SimpleGraph,
                                subgraph_species::Array{Symbol},
                                graph::SimpleGraph,
                                graph_species::Array{Symbol})::Dict{Int, Array{Bool, 2}}
        @debug "Building metadata dictionaries."
        𝒫s = DataFrame(index = 1:nv(subgraph), species = subgraph_species, degree = degree(subgraph))
        𝒫g = DataFrame(index = 1:nv(graph), species = graph_species, degree = degree(graph))
        @debug "Getting initial M and adjacency matrices."
        M° = correspondence_matrix(subgraph, 𝒫s, graph, 𝒫g)
        𝒜s = Array{Bool}(LinAlg.adjacency_matrix(subgraph))
        𝒜g = Array{Bool}(LinAlg.adjacency_matrix(graph))
        # Perform depth-first search
        ℳ = ullmann_DFS(M°, 𝒜s, 𝒜g)
        # Return solutions
        return ℳ
    end


## API


    @doc raw"""
        function subgraph_isomorphisms(search_moiety::Crystal,
            parent_structure::Crystal)::Union{Array{Array{Bool, 2}}, Nothing}
        Returns all bijections that map the bonding network from search_moiety
        onto parent_structure via Ullmann's algorithm.
    """
    function subgraph_isomorphisms(search_moiety::Crystal,
                                   parent_structure::Crystal)::Dict{Int, Array{Bool, 2}}
        @debug "Running Ullmann's algorithm to obtain bijections of $(search_moiety.name) into $(parent_structure.name)"
        return ullmann_bijections(search_moiety.bonds, search_moiety.atoms.species,
                                  parent_structure.bonds, parent_structure.atoms.species)
    end


    @doc raw"""
        function subgraph_find_replace(search_moiety::Crystal, replace_moiety::Crystal,
    							   parent_structure::Crystal)::Array{Array{Bool, 2}}
    	Returns a Crystal wherein all instances of search_moiety from parent_structure
    	are replaced by replace_moiety.  Uses Ullmann's algorithm for matching.
    """
    function subgraph_find_replace(search_moiety::Crystal, replace_moiety::Crystal,
    							   parent_structure::Crystal)::Array{Array{Bool, 2}}
    	matches = subgraph_isomorphisms(search_moiety, parent_structure)
    	return nothing
    end


    export subgraph_isomorphisms, subgraph_find_replace

end
