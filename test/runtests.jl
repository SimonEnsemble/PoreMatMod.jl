testfiles = [
    "moiety.jl"
    "ullmann.jl"
    "findreplace.jl"
]

using Test
using MOFun

@info "\n\n\t\tMOFun\n\n "

for testfile ∈ testfiles
    @info "Running test/$testfile"
    @time include(testfile)
end

@info "Done."
