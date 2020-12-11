TRAVIS=false # set false for better UX in local testing. must set true for CI!

testfiles = [
    "moiety.jl"
    "Ullmann.jl"
    "findreplace.jl"
]

## dependencies and logging
using Test
if !TRAVIS
    using Logging, Revise
    global_logger(ConsoleLogger(stdout, Logging.Info))
end

## banner
@info "\n\n\t\tMOFun Tests\n\n "

## race condition -> workers think modules not loaded. solution: just load again.
for _ in 1:2
    try
        using MOFun
    catch
    end
end

if !TRAVIS
    dirs = ["test", "data"]
else
    dirs = ["data"]
end
set_path_to_data(joinpath(pwd(), dirs...), relpaths=true, print=true)

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
if TRAVIS
    include.(testfiles)
else
    runtest.(testfiles)
end
