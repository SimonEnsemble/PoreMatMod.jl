# runtests.jl

## dependencies and logging
using Logging, Test, LightGraphs, LinearAlgebra, Xtals, Revise
global_logger(ConsoleLogger(stdout, Logging.Info))

## banner
@info "\n\n\t\tMOFun Tests\n\n "

## race condition -> workers think modules not loaded. solution: just load again.
for _ in 1:2
    try
        using PorousMaterials, MOFun
    catch
    end
end

## gives timing data for tests, plus some formatting
function runtest(testfile::String)
    @info "Testing $(testfile)"
    try
        @timev include(testfile)
    catch exception
        @error "Exception in $(testfile)" exception
    end
    println("")
end

## Run tests
runtest.([
    "moiety.jl"
    "Ullmann.jl"
    "findreplace.jl"
])
