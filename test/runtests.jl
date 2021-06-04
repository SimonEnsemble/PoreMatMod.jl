testfiles = [
    "moiety.jl"
    "ullmann.jl"
    "findreplace.jl"
]

using Test, Documenter

@info "\n\n\t\tMOFun\n\n "

using MOFun

doctest(MOFun, fix=true)

for testfile ∈ testfiles
    @info "Running test/$testfile"
    @time include(testfile)
end

@info "Done."
