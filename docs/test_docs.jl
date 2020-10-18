using Documenter
for _ in 1:2
    try
        using PorousMaterials, MOFun
    catch
    end
end

makedocs(
    root = joinpath(dirname(pathof(MOFun)), "..", "docs"),
    modules = [PorousMaterials],
    sitename = "MOFun.jl",
    clean = true,
    pages = [
        "MOFun" => "index.md",
        "How-To" => "guides/how2.md",
        "Manual" => [
            "Crystals and Moieties" => "manual/inputs.md",
            "Substructure Search" => "manual/find.md",
            "Substructure Find/Replace" => "manual/replace.md",
        ],
        "MOFunGO" => "mofungo.md"
    ],
    format = Documenter.HTML(assets = ["assets/flux.css"])
)
