using Documenter, PoreMatMod

makedocs( # to test docs, run this call to `makedocs` in the REPL
    root = joinpath(dirname(pathof(PoreMatMod)), "..", "docs"),
    modules = [PoreMatMod],
    sitename = "PoreMatMod.jl",
    clean = true,
    pages = [
        "PoreMatMod" => "index.md",
        "Getting Started" => "manual/start.md",
        "Loading Data" => "manual/inputs.md",
        "Substructure Search" => "manual/find.md",
        "Substructure Find/Replace" => "manual/replace.md",
        "PoreMatModGO" => "manual/PoreMatModGO.md",
        "Collaborate" => "collab.md",
        "Examples" => "manual/examples.md",
        "Hypothetical Structures" => "examples/make_hypothetical_MOF.jl.html",
        "Missing Atoms" => "examples/correct_missing_Hs.jl.html",
        "XRD Artifcats" => "examples/disorder_and_guests.jl.html",
        "Engineered Defects" => "examples/missing_linker_defect.jl.html"
    ],
    format = Documenter.HTML(assets = ["assets/flux.css"])
)

deploydocs(repo = "github.com/SimonEnsemble/PoreMatMod.jl.git")
