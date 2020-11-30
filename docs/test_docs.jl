using Documenter, Revise
for _ in 1:2
    try
        using MOFun
    catch
    end
end

makedocs(
    root = joinpath(dirname(pathof(MOFun)), "..", "docs"),
    modules = [MOFun],
    sitename = "MOFun.jl",
    clean = true,
    pages = [
        "MOFun" => "index.md",
        "Getting Started" => "manual/start.md",
        "Loading Data" => "manual/inputs.md",
        "Substructure Search" => "manual/find.md",
        "Substructure Find/Replace" => "manual/replace.md",
        "Examples" => "manual/examples.md",
        "MOFunGO" => "manual/mofungo.md",
        "Collaborate" => "collab.md"
    ],
    format = Documenter.HTML(assets = ["assets/flux.css"])
)
