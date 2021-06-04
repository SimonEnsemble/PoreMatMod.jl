testfiles = [
    "moiety.jl"
    "ullmann.jl"
    "findreplace.jl"
]

using Test, Documenter, FIGlet

FIGlet.render("PoreMatMod.jl", FIGlet.availablefonts()[35])

using MOFun

for testfile ∈ testfiles
    @info "Running test/$testfile"
    @time include(testfile)
end

@time doctest(MOFun)

@info "Done."
