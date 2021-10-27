module PoreMatMod

using StatsBase, Graphs, MetaGraphs, LinearAlgebra, DataFrames, Reexport
@reexport using Xtals

import Base.(∈), Base.show, Base.replace


function __init__()
    rc[:paths][:moieties] = joinpath(rc[:paths][:data], "moieties")
    rc[:r_tag] = '!'
    add_bonding_rules(tagged_bonding_rules())
end


export
    # search.jl
    Search, substructure_search, nb_isomorphisms, nb_locations, nb_ori_at_loc, isomorphic_substructures,

    # replace.jl
    substructure_replace, 

    # moiety.jl
    moiety

include("Ullmann.jl")
include("moiety.jl")
include("misc.jl")
include("search.jl")
include("replace.jl")
include("examples.jl")

end
