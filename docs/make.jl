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
        "Examples" => "manual/examples.md",
        "PoreMatModGO" => "manual/PoreMatModGO.md",
        "Collaborate" => "collab.md"
    ],
    format = Documenter.HTML(assets = ["assets/flux.css"])
)

deploydocs(repo = "github.com/SimonEnsemble/PoreMatMod.jl.git")
