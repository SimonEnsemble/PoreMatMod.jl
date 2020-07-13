#!/usr/bin/env julia
# ullmann_demo.jl
# Adrian Henle, 2020
# Demonstrates the function of MOFun's Ullmann bijection finder


using PorousMaterials, Logging
include("ullmann.jl")
include("moiety.jl")
global_logger(Logging.ConsoleLogger(stdout, Logging.Debug))


function ullmann_demo()
	@debug "Loading crystal and building bonds graph."
    xtal = moiety("DCM")
	strip_numbers_from_atom_labels!(xtal)
    #infer_bonds!(xtal, true)
	@debug "Crystal: $(xtal)"

	@debug "Loading search moiety."
    moty = moiety("methylene")
	@debug "Moiety: $(moty)"

    @debug "Finding subgraph isomorphisms."
    matches = Ullmann.subgraph_isomorphisms(moty, xtal)
	@debug "Results: $(matches)"

    # Replace with new Moiety
## TODO implement subgraph_find_replace
end


if(joinpath(pwd(),PROGRAM_FILE)==@__FILE__)||((@isdefined Atom)&&typeof(Atom)==Module)
	println()
	@debug "MAIN START"
	ullmann_demo()
	@debug "MAIN END"
end
