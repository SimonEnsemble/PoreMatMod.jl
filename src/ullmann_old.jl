# ullmann.jl
# Adrian Henle, 2020

@doc raw"""
    Ullmann algorithm code for PorousMaterials.jl
"""
module Ullmann
    export subgraph_isomorphisms, subgraph_find_replace
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
        # Allocate M at correct size (rows correspond to subgraph nodes, columns to graph nodes)
        M = zeros(Bool, nv(subgraph), nv(graph))
        for i ∈ 1:nv(subgraph) # Loop over rows
            for j ∈ 1:nv(graph) # Loop over columns
@debug "($i,$j) $(𝒫g.degree[j]) ≥ $(𝒫s.degree[i]) && $(𝒫g.species[j]) == $(𝒫s.species[i]) : $(𝒫g.degree[j] ≥ 𝒫s.degree[i] && 𝒫g.species[j] == 𝒫s.species[i]))"
                # Record Bool for each (i,j): true if atom species match and graph node degree is sufficient.
                M[i, j] = 𝒫g.degree[j] ≥ 𝒫s.degree[i] && 𝒫g.species[j] == 𝒫s.species[i]
            end
        end
        @debug "M: $(M)"
        return M
    end


    @doc raw"""
        function suppose_correspondence(M0::Array{Bool, 2}, subgraph_node::Int,
                                        graph_node::Int)::Array{Bool, 2}
        Generates a new matrix from M0 by supposing x corresponds to y
    """
    function suppose_correspondence(M0::Array{Bool, 2}, subgraph_node::Int,
                                    graph_node::Int)::Array{Bool, 2}
        @debug "Making supposition S[$(subgraph_node)]↔G[$(graph_node)]"
        M = deepcopy(M0) # Do not mutate parent M
        M[subgraph_node, :] .= false # zero out row of subgraph node
        M[:, graph_node] .= false # zero out column of graph node
        M[subgraph_node, graph_node] = true # assign correspondence at intersection
        return M
    end


    @doc raw"""
        function candidate_list(M::Array{Bool,2},node_index::Int)::Array{Int}
        Gives the list of candidates for mapping a subgraph node according to M.
    """
    function candidate_list(M::Array{Bool,2}, node_index::Int)::Array{Int}
        @debug "Getting candidate list for subgraph node $(node_index)"
        # Get Bools for correspondence w/ subgraph node
        logicals = M[node_index, :]
        @debug "Logicals: $(logicals)"
        # Translate to list of graph node indices
        return [index for (index, logical) ∈ enumerate(logicals) if logical]
    end


    @doc raw"""
        function validate_M(M::Array{Bool, 2})::Bool
        Determines if M has no empty candidate lists for subgraph nodes.
    """
    function validate_M(M::Array{Bool, 2})::Bool
        @debug "Validating M: $(M)"
        for S_node ∈ 1:size(M, 1) # loop over subgraph nodes
            if length(candidate_list(M, S_node)) == 0
                return false # if any subgraph node cannot be assigned, M is the intersection of zero solutions.
            end
        end
        return true # M may be the intersection of one or more solutions.
    end


    @doc raw"""
        function is_solution(M::Array{Bool, 2})::Bool
        Determines if M is a 1-to-1 mapping solution
    """
    function is_solution(M::Array{Bool, 2})::Bool
        @debug "Testing for solution: $(M)"
        for i ∈ 1:size(M, 1)
            if length(candidate_list(M, i)) ≠ 1
                @debug "Not a solution."
                return false # Only a solution if every subgraph node has exactly one possible correspondence
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
        # Return the index of each node sharing an edge with node w
        return [v for v ∈ 1:size(A, 1) if A[v, w]]
    end


    @doc raw"""
        function refine_M!(As::Array{Bool, 2}, Ag::Array{Bool, 2}, M::Array{Bool, 2})
        Performs Ullman algorithm refinement until M is stable
    """
    function refine_M!(As::Array{Bool, 2}, Ag::Array{Bool, 2}, M::Array{Bool, 2})
        @debug "Refining M: $(M)"
        for i ∈ 1:size(As, 1)*size(Ag, 1)+1 # Max num of iterations: product of graph node counts + 1
            M_altered = false # Tracks if M has been altered in the current iteration
            for y ∈ size(M, 1) # Loop over subgraph nodes
                for x ∈ candidate_list(M, y) # Loop over candidates of subgraph node
                    for z ∈ neighbors(y, As) # Loop over neighbors of subgraph node
                        if length(intersect(candidate_list(M, z), neighbors(x, Ag))) == 0
                            # M[y,x] can only be true if at least one neighbor of y∈subgraph has possible correspondence to at least one neighbor of x∈graph
                            M[y, x] = false
                            @debug "M altered at [$(y), $(x)]."
                            M_altered = true
                        end
                    end
                end
            end
            if !M_altered # When M not altered by one full refinement cycle, it is stable (refinement done)
                return M
            end
        end
    end


    @doc raw"""
        function ullmann_DFS(M::Array{Bool, 2}, AS::Array{Bool, 2}, AG::Array{Bool, 2},
                       i°::Int=1, j°::Int=1)::Union{Array{Array{Bool, 2}}, Nothing}
        Performs depth-first search for Ullmann's algorithm.
    """
    function ullmann_DFS(M::Array{Bool, 2}, AS::Array{Bool, 2}, AG::Array{Bool, 2},
                       i°::Int=1, dfs_rec_lvl::Int=-1)::Dict{Int, Array{Bool, 2}}
        dfs_rec_lvl += 1 # tracker for depth of search
        @debug "DFS Recursion Level: $dfs_rec_lvl"
        ℳ = Dict{Int, Array{Bool, 2}}() # solution set
        if validate_M(M) # Only search if M is intersection of ≥1 possible solutions
            for i ∈ i°:size(M, 1) # Looping over subgraph nodes (starting from lowest uncertain index)
                @debug "i: $i"
                if length(candidate_list(M, i)) > 1 # Only work on this subgraph node if it has multiple possible correspondences
## TODO refactor (new) candidate_list for rapid return here
                    for j ∈ 1:size(M, 2) # Looping over graph nodes
                        @debug "j: $j"
                        if M[i, j] # Skip over non-possible correspondences
                            M′ = suppose_correspondence(M, i, j) # Eliminate all other correspondences for the (i,j) pair
                            refine_M!(AS, AG, M′) # Eliminate possibilities incompatible with supposition
                            if is_solution(M′) # If solution reached, record and test other correspondences for current i
                                @debug "Appending solution to ℳ: $(M′)"
                                ℳ[length(ℳ)+1] = copy(M′)
                                @debug "Solutions:" ℳ
                            elseif validate_M(M′) # If not solution but may lead to one, search recursively
                                @debug "Continuing search..."
                                merge(ℳ, ullmann_DFS(M′, AS, AG, i+1, dfs_rec_lvl)) # THIS ISN'T GOOD
                                # merge() can overwrite solutions!
## TODO refactor to Array{Array{Bool}}
                                dfs_rec_lvl -= 1
## TODO don't need to track explicitly (depth = i)
                            else
                                # Reached a state where M′ is not a solution and does not lead to one.
                                # Continue checking other values of j
                                @debug "Reached dead end."
                            end
                        end
                    end
                end
            end
        end
        dfs_rec_lvl -= 1
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
        # Dicts for holding graph metadata: atomic species, node sort keys
        𝒫s = DataFrame(index = 1:nv(subgraph), species = subgraph_species, degree = degree(subgraph))
        𝒫g = DataFrame(index = 1:nv(graph), species = graph_species, degree = degree(graph))
        @debug "Getting initial M and adjacency matrices."
        # M° encodes intersection of every possible solution
        M° = correspondence_matrix(subgraph, 𝒫s, graph, 𝒫g)
        # Symmetric square matrices encoding graph edges
        𝒜s = Array{Bool}(LinAlg.adjacency_matrix(subgraph))
        𝒜g = Array{Bool}(LinAlg.adjacency_matrix(graph))
        @debug "Starting depth-first search."
        # Return depth-first search results
        return ullmann_DFS(M°, 𝒜s, 𝒜g)
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
        # Defensive input validation
        @assert ne(search_moiety.bonds) > 0 "Search moiety must have bonds defined!"
        @assert ne(parent_structure.bonds) > 0 "Parent structure must have bonds defined!"
        @assert nv(search_moiety.bonds) ≤ nv(parent_structure.bonds) "Search moiety must be smaller than parent structure!"
        # Wrapping of search method
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
        # TO BE IMPLEMENTED after subgraph_isomorphisms() is working
    	matches = subgraph_isomorphisms(search_moiety, parent_structure)
    	return nothing
    end
end
