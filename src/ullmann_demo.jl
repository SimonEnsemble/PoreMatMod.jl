#!/usr/bin/env julia
# ullmann_demo.jl
# Adrian Henle, 2020
# Demonstrates the function of MOFun's Ullmann bijection finder

## TODO
"""
    # Implement find_replace()
"""
## #

using PorousMaterials
include("ullmann.jl")
include("moiety.jl")


function ullmann_demo()
    # Build parent Crystal
    Co2m_dobdc = Crystal("KOSKIO_clean.cif")
    infer_bonds!(Co2m_dobdc, true)

    # Build search Moiety
    m_dobdc = moiety("data/fragments/m-dobdc")

    # Search for search Moiety in parent Crystal
    matches = Ullmann.subgraph_isomorphisms(m_dobdc, Co2m_dobdc)

    # Replace with new Moiety
    f2_m_dobdc = moiety("data/fragments/F2-m-dobdc")
## TODO implement subgraph_find_replace
end

ullmann_demo()

"""if (abspath(PROGRAM_FILE) == @__FILE__) || ((@isdefined Atom) && typeof(Atom) == Module)
    ullmann_demo()
end"""
