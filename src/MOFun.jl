module MOFun

using StatsBase, LightGraphs, MetaGraphs, LinearAlgebra, DataFrames, Reexport
@reexport using Xtals

import Base.(∈), Base.show


function __init__()
    rc[:paths][:moieties] = joinpath(rc[:paths][:data], "moieties")
    rc[:r_tag] = '!'
    add_bonding_rules(tagged_bonding_rules())
end


export
    # Xtals
    rc,

    # findreplace.jl
    substructure_search, SearchResult, Query, Search, nb_isomorphisms,
    nb_locations, nb_configs_at_loc, replace,

    # moiety.jl
    moiety

include("Ullmann.jl")
include("moiety.jl")
include("findreplace.jl")

end
