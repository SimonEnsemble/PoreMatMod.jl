module PoreMatMod

using StatsBase, Graphs, MetaGraphs, LinearAlgebra, DataFrames, Pluto, Reexport
@reexport using Xtals

import Base.(∈), Base.show, Base.replace


function __init__()
    add_bonding_rules(tagged_bonding_rules())
end


"""
    PoreMatModGO()

Launches the GUI Pluto notebook.
"""
function PoreMatModGO()
    # notebook path
    pmmg_ntbk = joinpath(pathof(PoreMatMod), "..", "PoreMatModGO.jl")
    # run the notebook in Pluto
    Pluto.run(notebook=pmmg_ntbk)
end


export
    # search.jl
    Search, substructure_search, nb_isomorphisms, nb_locations, nb_ori_at_loc, isomorphic_substructures,

    # replace.jl
    substructure_replace, 

    # moiety.jl
    moiety,

    # misc.jl
    PoreMatModGO

include("Ullmann.jl")
include("moiety.jl")
include("search.jl")
include("replace.jl")
include("ExampleHelpers.jl")

end
