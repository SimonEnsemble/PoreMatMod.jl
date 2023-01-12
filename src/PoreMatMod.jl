module PoreMatMod

using DataFrames, Graphs, LinearAlgebra, MetaGraphs, Pluto, Reexport, StatsBase
@reexport using Xtals
using PrecompileSignatures: @precompile_signatures

import Base.(∈), Base.show, Base.replace

__init__() = add_bonding_rules(tagged_bonding_rules())

export
    # search.jl
    Search,
    substructure_search,
    nb_isomorphisms,
    nb_locations,
    nb_ori_at_loc,
    isomorphic_substructures,

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
include("misc.jl")

@precompile_signatures(PoreMatMod)

end
