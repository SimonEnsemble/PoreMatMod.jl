#!/usr/bin/env julia
# ullmann_demo.jl
# Adrian Henle, 2020
# Demonstrates the function of MOFun's Ullmann bijection finder

## TODO
"""
    # Refactor ullmann_bijections() call as xyz_find()
    # Implement xyz_find() API in MOFun.jl
    # Implement xyz_find_replace()
"""
## #

#using MOFun
include("MOFun.jl")


function ullmann_demo()
    # Build parent Crystal
    Co2m_dobdc = Crystal("KOSKIO_clean.cif")
    infer_bonds!(Co2m_dobdc, true)

    # Build search Moiety
    m_dobdc = moiety("data/fragments/m-dobdc.xyz")

    # Search for search Moiety in parent Crystal
    matches = subgraph_isomorphisms(m_dobdc, Co2m_dobdc)

    # Replace with new Moiety
    f2_m_dobdc = moiety("data/fragments/F2-m-dobdc.xyz")
## TODO implement subgraph_find_replace
end

if (abspath(PROGRAM_FILE) == @__FILE__) || ((@isdefined Atom) && typeof(Atom) == Module)
    using Profile
    @profile ullmann_demo()
end
