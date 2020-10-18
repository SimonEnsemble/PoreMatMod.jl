using Documenter
using PorousMaterials, MOFun

makedocs(
    root = joinpath(dirname(pathof(MOFun)), "..", "docs"),
    modules = [PorousMaterials],
    sitename = "MOFun.jl",
    clean = true,
    pages = [
            "MOFun" => "index.md",
            "How-To" => "guides/how2.md",
            "Manual" => [
                "Crystals and Moieties" =>
                    "manual/inputs.md",

                "Substructure Search" =>
                    "manual/find.md",

                "Substructure Find/Replace" =>
                    "manual/replace.md",
                ],
            "MOFunGO" => "mofungo.md"
            ],
    format = Documenter.HTML(assets = ["assets/flux.css"])
)

deploydocs(
    repo = "github.com/SimonEnsemble/PorousMaterials.jl.git"
 #     push_preview=false,
 #     deps = Deps.pip("mkdocs", "mkdocs-material", "pymdown-extensions") # These are dependencies for the site, not the package
)
