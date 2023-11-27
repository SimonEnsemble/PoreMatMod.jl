module MarkdownExt


using PoreMatMod

using Graphs, Markdown, Reexport, MetaGraphs, LinearAlgebra, Xtals

include("shared/required_files.jl")

__init__() = check_example_data()

function input_file_message()
    return md"""
!!! note \"input files for the example Pluto notebooks\"
    if the input files required for the example Pluto notebooks are not present in the correct folders, `ExampleHelpers` automatically copies the required input files from the `examples/data` directory of the `PoreMatMod.jl` source code to the folders `rc[:paths][:crystals]` and `rc[:paths][:moieties]`. all input files for the examples are also on Github [here](https://github.com/SimonEnsemble/PoreMatMod.jl/tree/master/examples/data).

    n.b. you may change the folders from which `PoreMatMod.jl` reads input files by setting `rc[:paths][:crystals]` and `rc[:paths][:moieties]` as the desired path. for example, if you desire to store your crystal structures in a folder `~/my_xtals/` (a folder in your home directory), set:
    ```julia
    rc[:paths][:crystals] = joinpath(homedir(), \"my_xtals\").
    ```
"""
end

function xtal_folder_message()
    return md"""
📕 folder from which `PoreMatMod.jl` reads `.cif` files that represent crystal structures:
"""
end

function moiety_folder_message()
    return md"""
📕 folder from which `PoreMatMod.jl` reads `.xyz` files that represent fragments/moities:
"""
end

function fragment_construction_note()
    return md"""
!!! note \"how can we construct the query/replacement fragments?\"
    two options we use: (1) use Avogadro as a molecule builder/editor and export it as `.xyz` or (2) cut the appropriate fragment out of the MOF crystal structure in the `.cif` file using e.g. iRASPA.
"""
end

write_cif_message() = md"""
           write the child crystal structure to file for downstream molecular simulations
           """

# function to print the contents of a moiety file
function display_query_or_replacement_file(filename::String)
    filename = joinpath(rc[:paths][:moieties], filename)

    println("contents of: ", filename, "\n")
    open(filename, "r") do io
        return print(read(io, String))
    end
end

export display_query_or_replacement_file,
    write_cif_message,
    xtal_folder_message,
    moiety_folder_message,
    fragment_construction_note,
    input_file_message,
    check_example_data

end
