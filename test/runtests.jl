testfiles = [
    "moiety.jl"
    "ullmann.jl"
    "findreplace.jl"
    "examples.jl"
]

@assert VERSION.major == 1
@assert VERSION.minor ≥ 6

using Test, Documenter, PoreMatMod
PoreMatMod.banner()

for testfile ∈ testfiles
    @info "Running test/$testfile"
    @time include(testfile)
end

if VERSION.minor ≥ 7
    @time doctest(PoreMatMod)
end

@info "Done."
