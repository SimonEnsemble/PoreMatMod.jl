module MOFun

using LightGraphs
using LinearAlgebra
#using Printf
using DataFrames

using Xtals
export Crystal

include("Ullmann.jl")
include("moiety.jl")
include("findreplace.jl")


"""
	set_path_to_moieties("../moieties")

Assigns the path for locating .xyz moiety input files.

# Arguments
- `path::String` the location of the moiety input files
- `print::Bool` set true to print the `Xtals` and `MOFun` data path variables
"""
function set_path_to_moieties(path::String; print::Bool=false)
	global PATH_TO_MOIETIES = path
	if print
		print_file_paths()
	end
end

"""
	set_path_to_data("../data")

Assigns the path for locating general data.

# Arguments
- `path::String` the location of the general data folder
- `relpath_xtals::Bool` set true to automatically assign `Xtals.PATH_TO_CRYSTALS` relative to `Xtals.PATH_TO_DATA`
- `relpath_motys::Bool` set true to automatically asign `MOFun.PATH_TO_MOIETIES` relative to `Xtals.PATH_TO_DATA`
- `print::Bool` set true to print the `Xtals` and `MOFun` data path variables
"""
function set_path_to_data(path::String; relpath_xtals::Bool=false, relpath_motys::Bool=false, print::Bool=false)
	Xtals.set_path_to_data(path, relpath_xtals)
	if relpath_motys
		set_path_to_moieties(joinpath(PATH_TO_DATA, "moieties"))
	end
	if print
		print_file_paths()
	end
end

"""
	print_file_paths()

Prints the `Xtals` and `MOFun` data path variables.
"""
function print_file_paths()
	Xtals.print_file_paths()
	println("moieties (.xyz): ", PATH_TO_MOIETIES)
end








end
