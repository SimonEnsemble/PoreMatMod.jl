using Revise, Logging, Test
global_logger(ConsoleLogger(stdout, Logging.Info))
using LightGraphs, LinearAlgebra, Printf
using PorousMaterials, MOFun, Ullmann, Moiety
function runtest(testfile::String)
    @info "Testing $(testfile)"
    try
        @timev include(testfile)
    catch exception
        @error "Exception in $(testfile)" exception
    end
end

testfiles = [#"alignment_operations.jl"
             #"ring_constructor.jl"
             "Moiety.jl"
             "Ullmann.jl"
             "MOFun.jl"
             ]

[runtest(testfile) for testfile in testfiles]
