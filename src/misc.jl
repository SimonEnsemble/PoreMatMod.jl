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


"""
    banner()

Displays the ASCII banner for the package using FIGlet
"""
banner() = FIGlet.render("PoreMatMod.jl", FIGlet.availablefonts()[35])