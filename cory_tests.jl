using PorousMaterials, LightGraphs

include("src/moiety.jl")
fragment = moiety("p-phenylene")

xtal = Crystal("IRMOF-1.cif")
strip_numbers_from_atom_labels!(xtal)
infer_bonds!(xtal, false)
xtal.bonds

subgraph = fragment.bonds
graph = xtal.bonds

push!(LOAD_PATH, joinpath(pwd(), "src"))
using Ullmann

solns = find_subgraph_isomorphisms(fragment.bonds, fragment.atoms.species, xtal.bonds, xtal.atoms.species)
println("found solns: ", solns)
