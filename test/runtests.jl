testfiles = [
    "moiety.jl"
    "ullmann.jl"
    "findreplace.jl"
    "examples.jl"
]

using Test, Documenter, FIGlet

FIGlet.render("PoreMatMod.jl", FIGlet.availablefonts()[35])

using PoreMatMod

for testfile ∈ testfiles
    @info "Running test/$testfile"
    @time include(testfile)
end

@time doctest(PoreMatMod)

@info "Done."
