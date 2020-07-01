#!/usr/bin/env julia
# ullmann_demo.jl
# Adrian Henle, 2020
# Demonstrates the function of MOFun's Ullmann bijection finder

## TODO
"""
    # Implement find_replace()
"""
## #

using PorousMaterials, Logging
include("ullmann.jl")
include("moiety.jl")
global_logger(Logging.ConsoleLogger(stdout, Logging.Debug))


function ullmann_demo()
	@debug "Loading crystal and building bonds graph."
    Co2m_dobdc = Crystal("KOSKIO_clean.cif")
    infer_bonds!(Co2m_dobdc, true)
	@debug "Crystal: $(Co2m_dobdc)"

	@debug "Loading search moiety."
    m_dobdc = moiety("data/moieties/m-dobdc")
	@debug "Moiety: $(m_dobdc)"

    @debug "Finding subgraph isomorphisms."
    matches = Ullmann.subgraph_isomorphisms(m_dobdc, Co2m_dobdc)
	@debug "Results: $(matches)"

    # Replace with new Moiety
    f2_m_dobdc = moiety("data/moieties/F2-m-dobdc")
## TODO implement subgraph_find_replace
end

if(joinpath(pwd(),PROGRAM_FILE)==@__FILE__)||((@isdefined Atom)&&typeof(Atom)==Module)
	println()
	@debug "MAIN START"
	ullmann_demo()
	@debug "MAIN END"
end
