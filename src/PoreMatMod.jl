module PoreMatMod

using DataFrames, Graphs, LinearAlgebra, MetaGraphs, Reexport, StatsBase
@reexport using Xtals
# using PrecompileSignatures: @precompile_signatures

import Base.(∈), Base.show, Base.replace

const BANNER = String(read(joinpath(dirname(pathof(PoreMatMod)), "banner.txt")))
banner() = println(BANNER)

function __init__()
    rc[:r_tag] = '!'
    rc[:paths][:moieties] = joinpath(rc[:paths][:data], "moieties")
    add_bonding_rules(tagged_bonding_rules())
    return
end

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

    # MarkdownExt
    input_file_message

include("Ullmann.jl")
include("moiety.jl")
include("search.jl")
include("replace.jl")

# @precompile_signatures(PoreMatMod)

end
