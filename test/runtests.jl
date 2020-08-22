using Logging, Test, LightGraphs, LinearAlgebra, Revise
global_logger(ConsoleLogger(stdout, Logging.Info))
using PorousMaterials, MOFun, Ullmann

function runtest(testfile::String)
    @info "Testing $(testfile)"
    try
        @timev include(testfile)
    catch exception
        @error "Exception in $(testfile)" exception
    end
    println("")
end

testfiles = [#"alignment_operations.jl"
             #"ring_constructor.jl"
             "moiety.jl"
             "Ullmann.jl"
             "MOFun.jl"
             ]

@info "\n\n\t\tMOFun Tests\n\n "
[runtest(testfile) for testfile in testfiles]
