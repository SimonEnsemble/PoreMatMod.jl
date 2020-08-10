using Logging
global_logger(ConsoleLogger(stdout, Logging.Info))
using PorousMaterials, LightGraphs, MOFun, Ullmann, Test, LinearAlgebra, Printf
function runtest(testfile::String)
    @info "Testing $(testfile)"
    try
        @timev include(testfile)
    catch exception
        @error "Exception in $(testfile)" exception
    end
end

testfiles = [#"alignment_operations.jl"
             #"ring_constructor"
             "Ullmann.jl"
             "MOFun.jl"
             ]

[runtest(testfile) for testfile in testfiles]
