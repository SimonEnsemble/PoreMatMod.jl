module PlutoExt

using PoreMatMod

"""
    PoreMatModGO()

Launches the GUI Pluto notebook.
"""
function PoreMatModGO()
    try
        # notebook path
        pmmg_ntbk = joinpath(pathof(PoreMatMod), "..", "PoreMatModGO.jl")
        # run the notebook in Pluto
        Pluto.run(; notebook=pmmg_ntbk)
    catch
        # download as a temporary file in case of access issues
        pmmg_ntbk = download(
            "https://raw.githubusercontent.com/SimonEnsemble/PoreMatMod.jl/master/src/PoreMatModGO.jl"
        )
        Pluto.run(; notebook=pmmg_ntbk)
    end
end

export PoreMatModGO

end
