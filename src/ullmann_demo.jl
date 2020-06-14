# ullmann_demo.jl
# Adrian Henle, 2020
# Demonstrates the function of MOFun's Ullmann bijection finder

if abspath(PROGRAM_FILE) == @__FILE__
    # Build parent Crystal
    Co2m_dobdc = Crystal("KOSKIO_clean.cif")
    infer_bonds!(Co2m_dobdc, true)

    # Build search Moiety
    m_dobdc = moiety("data/fragments/m-dobdc.xyz")

    # Search for search Moiety in parent Crystal
    matches = ullmann_bijections(m_dobdc, Co2m_dobdc)

    # Replace with new Moiety
    f2_m_dobdc = moiety("data/fragments/F2-m-dobdc.xyz")
    ## TODO

    exit(0)
end
