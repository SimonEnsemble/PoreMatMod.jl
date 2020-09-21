# runtests.jl

## dependencies and logging
using Logging, Test, LightGraphs, LinearAlgebra, PorousMaterials, Revise
global_logger(ConsoleLogger(stdout, Logging.Info))

## banner
@info "\n\n\t\tMOFun Tests\n\n "

## silly hack to handle dumb errors loading local modules
try
    using Ullmann
catch
    using Ullmann
end
try
    using MOFun
catch
    using MOFun
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

## Run the tests!
runtest.([
    "moiety.jl"
    "Ullmann.jl"
    "MOFun.jl"
])
