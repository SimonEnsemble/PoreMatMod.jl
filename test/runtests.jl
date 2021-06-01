testfiles = [
    "moiety.jl"
    "ullmann.jl"
    "findreplace.jl"
]

## dependencies and logging
using Test

## banner
@info "\n\n\t\tMOFun\n\n "

using MOFun

for testfile ∈ testfiles
    @info "Running test/$testfile"
    @time include(testfile)
end

@info "Done."
