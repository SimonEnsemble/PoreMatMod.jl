#!/usr/bin/env julia
# ullmann_demo.jl
# Adrian Henle, 2020
# Demonstrates the function of MOFun's Ullmann bijection finder

function ullmann_demo()
    # Build parent Crystal
    Co2m_dobdc = Crystal("KOSKIO_clean.cif")
    infer_bonds!(Co2m_dobdc, true)

    # Build search Moiety
    m_dobdc = moiety("data/fragments/m-dobdc.xyz")

    # Search for search Moiety in parent Crystal
    matches = ullmann_bijections(m_dobdc, Co2m_dobdc)
## TODO refactor as xyz_find

    # Replace with new Moiety
    f2_m_dobdc = moiety("data/fragments/F2-m-dobdc.xyz")
## TODO implement xyz_find_replace
end

if (abspath(PROGRAM_FILE) == @__FILE__) || ((@isdefined Atom) && typeof(Atom) == Module)

    ullmann_demo()

    exit(0)
end
