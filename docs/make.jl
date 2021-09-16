using Documenter, PoreMatMod, PlutoSliderServer

# run the Pluto example notebooks and export them to HTML for inclusion in the docs
cd("examples") # next line fails if not actually in examples/
PlutoSliderServer.export_directory() # runs each notebook in examples/ and exports to HTML
rm("index.html") # previous line makes additional table-of-contents page (not needed)
for file in readdir() # loop over files in examples/
    if length(split(file, ".html")) > 1 # select only the HTML files
        mv(file, "../docs/src/examples/$file", force=true) # move HTML files to docs source folder for build
    end
end
cd("..") # return to root for running makedocs

# build the docs
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
    format = Documenter.HTML(assets = ["assets/flux.css"]),
    doctest = false # doctests are run in testing; running them here is redundant and slow
)

# deploy the docs
deploydocs(repo = "github.com/SimonEnsemble/PoreMatMod.jl.git", push_preview=true)
