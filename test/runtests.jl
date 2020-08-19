using Revise, Logging, Test
global_logger(ConsoleLogger(stdout, Logging.Info))
#using LightGraphs, LinearAlgebra
using PorousMaterials, MOFun, Ullmann, Moiety
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
             "Moiety.jl"
             "Ullmann.jl"
             "MOFun.jl"
             ]

@info "\n\n\t\tMOFun Tests\n\n "
[runtest(testfile) for testfile in testfiles]
