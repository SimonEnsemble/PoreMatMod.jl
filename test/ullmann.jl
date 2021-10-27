module Ullmann_Test

using Test, Graphs, MetaGraphs, PoreMatMod


# add a list of edges to a graph
function add_edges!(g, edges)
    for edge in edges
        add_edge!(g, edge[1], edge[2])
    end
end


# make a graph with the listed labels and edges
function build_graph(labels, edges)
    g = MetaGraph()
    add_vertices!(g, length(labels))
    add_edges!(g, edges)
    return g, labels
end


@testset "Ullmann Tests" begin
    # (graph, labels)
    (g, lg) = build_graph([:A, :B, :B, :C, :A, :A], [(1,2),(2,3),(3,4),(4,2),(3,5),(5,6)])
    (s1, l1) = build_graph([:B, :C], [(1,2)])
    (s2, l2) = build_graph([:A, :B], [(1,2)])
    (s3, l3) = build_graph([:B, :B, :C], [(1,2),(2,3),(3,1)])
    (s4, l4) = build_graph([:B, :D, :C], [(1,2),(2,3),(3,1)])
    (s5, l5) = build_graph([:B, :A, :B, :C, :D], [(1,2),(1,3),(1,4),(1,5)])
    (s6, l6) = build_graph([:A, :A], [])
    (s7, l7) = build_graph([:A, :A, :B], [(1,2),(2,3)])

    @test PoreMatMod.find_subgraph_isomorphisms(g, lg, g, lg) == [[1, 2, 3, 4, 5, 6]]
    @test PoreMatMod.find_subgraph_isomorphisms(s1, l1, g, lg) == [[2, 4], [3, 4]]
    @test PoreMatMod.find_subgraph_isomorphisms(s2, l2, g, lg)== [[1, 2], [5, 3]]
    @test PoreMatMod.find_subgraph_isomorphisms(s3, l3, g, lg) == [[2, 3, 4], [3, 2, 4]]
    @test PoreMatMod.find_subgraph_isomorphisms(s4, l4, g, lg) == []
    @test PoreMatMod.find_subgraph_isomorphisms(s5, l5, g, lg) == []
    @test PoreMatMod.find_subgraph_isomorphisms(s7, l7, g, lg) == [[6, 5, 3]]
end
end
