# PMUllmann.jl
# Adrian Henle, 2020

"""
    Ullmann algorithm code for PorousMaterials.jl
"""
## TODO
"""
    # sort M° on node degree in ullmann_bijections()
    # validate on sufficient unique candidates?
    # other accelerations/boosts?
    # ullmann tests
    # API → MOFun.jl
"""
## #

module Ullmann

    using PorousMaterials, LightGraphs, DataFrames, LightGraphs.LinAlg, Logging
    global_logger(Logging.ConsoleLogger(stdout, Logging.Debug))

    @doc raw"""
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
        Gives the list of candidates for mapping the node according to M
        Argument `graph` must be a symbol, either :graph or :subgraph,
        to indicate to which graph the node belongs.
    """
        @debug "Getting candidate list for node $(node_index) of $(graph)"
        if graph == :subgraph
            logicals = M[node_index, :]
        elseif graph == :graph
            logicals = M[:, node_index]
        else
            @error "Must provide either :graph or :subgraph."
        end
        @debug "Logicals: $(logicals)"
        for (index, logical) ∈ enumerate(logicals)
            if logical
                @debug "$index, $logical"
            end
        end
        return [index for (index, logical) ∈ enumerate(logicals) if logical]
    end


    @doc raw"""
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
        Determines if M is a 1-to-1 mapping solution
    """
    function is_solution(M::Array{Bool, 2})::Bool
        @debug "Testing for solution: $(M)"
        for i ∈ 1:size(M, 1)
            if length(candidate_list(M, i, :subgraph)) ≠ 1
                return false
            end
        end
        for i ∈ 1:size(M, 2)
            if length(candidate_list(M, i, :graph)) > 1
                return false
            end
        end
        return true
    end


    @doc raw"""
        Returns list of neighbors of input node from its graph's adjacency matrix.
    """
    function neighbors(w::Int, A::Array{Bool, 2})::Array{Int}
        @debug "Finding neighbors of node $(w)."
        return [v for v ∈ 1:size(A, 1) if A[v, w]]
    end


    @doc raw"""
        Performs Ullman algorithm refinement until M is stable
    """
    function refine_M!(As::Array{Bool, 2}, Ag::Array{Bool, 2}, M::Array{Bool, 2})
        @debug "Refining M: $(M)"
        M_altered = false
        for y ∈ size(M, 1)
            for x ∈ candidate_list(M, y, :graph)
                for z ∈ neighbors(y, As)
                    if candidate_list(M, z, :subgraph)
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
        Performs depth-first search for Ullmann's algorithm.
    """
    function ullmann_DFS(M::Array{Bool, 2}, AS::Array{Bool, 2},
                         AG::Array{Bool, 2})::Union{Array{Array{Bool, 2}}, Nothing}
        ℳ = []
        if validate_M(M)
            for i ∈ 1:size(M, 1)
                if length(candidate_list(M, i, :graph)) > 1
                    for j ∈ 1:size(M, 2)
                        if M[i, j]
                            M′ = suppose_correspondence(M, i, j)
                            refine_M!(M′)
                            if is_solution(M′)
                                append!(ℳ, M′)
                            elseif validate_M(M′)
                                append!(ℳ, ullmann_DFS(M′))
                            else
                                return ℳ
                            end
                        end
                    end
                end
            end
            return ℳ
        end
    end


    @doc raw"""
        Returns all correspondence matrix bijections of subgraph to a subset of graph.
    """
    function ullmann_bijections(subgraph::SimpleGraph,
                                subgraph_species::Array{Symbol},
                                graph::SimpleGraph,
                                graph_species::Array{Symbol})::Union{Array{Array{Bool,2}}, Nothing}
        @debug "Building metadata dictionaries."
        𝒫s = DataFrame(index = 1:nv(subgraph), species = subgraph_species, degree = degree(subgraph))
        𝒫g = DataFrame(index = 1:nv(graph), species = graph_species, degree = degree(graph))
        @debug "Getting initial M and adjacency matrices."
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




    @doc raw"""
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


    @doc raw"""
    	Returns all bijections that map the bonding network from search_moiety onto
    	parent_structure via Ullmann's algorithm
    """
    function subgraph_isomorphisms(search_moiety::Crystal,
    							   parent_structure::Crystal)::Union{Array{Array{Bool, 2}}, Nothing}
        @debug "Running Ullmann's algorithm to obtain bijections of $(search_moiety.name) into $(parent_structure.name)"
        return ullmann_bijections(search_moiety.bonds, search_moiety.atoms.species,
    							  parent_structure.bonds, parent_structure.atoms.species)
    end

export subgraph_isomorphisms

end
