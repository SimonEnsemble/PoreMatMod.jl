using Logging

# list of required files for examples
required_files = Dict(:crystals => ["IRMOF-1.cif", "EMEHUB_C2H2.cif", "IRMOF-1_noH.cif", "UiO-66.cif"],
                      :moieties => ["2-!-p-phenylene.xyz", "2-acetylamido-p-phenylene.xyz", "1,4-C-phenylene_noH.xyz", "1,4-C-phenylene.xyz", "4-pyridyl.xyz",
					  "acetylene.xyz", "BDC.xyz", "disordered_ligand!.xyz", "formate_caps.xyz"]
                     )

# make sure directories are present and the right files for the examples
for file_type in [:moieties, :crystals]
    # make sure directories exist
    if ! isdir(rc[:paths][file_type])
        @warn "$(rc[:paths][file_type]) directory not present; creating it."
        mkpath(rc[:paths][file_type])
    end
    for required_file in required_files[file_type]
        where_it_shld_be = joinpath(rc[:paths][file_type], required_file)
        if ! isfile(where_it_shld_be)
            @warn "$where_it_shld_be not present; copying it from src."
            where_it_is = normpath(joinpath(pathof(PoreMatMod), "..", "..", "examples", "data", String(file_type), required_file))
            cp(where_it_is, where_it_shld_be)
        end
    end
end

function input_file_message() 
	return md"""
!!! note \"input files for the example Pluto notebooks\"
	if the input files required for the example Pluto notebooks are not present in the correct folders, the script `ExampleHelper.jl` automatically copies the required input files from the `examples/data` directory of the `PoreMatMod.jl` source code to the folders `rc[:paths][:crystals]` and `rc[:paths][:moieties]`. all input files for the examples are also on Github [here](https://github.com/SimonEnsemble/PoreMatMod.jl/tree/master/examples/data).

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

function write_cif_message()
	return md"""
write the child crystal structure to file for downstream molecular simulations
"""
end

# function to visualize a crystal in the notebook
function view_structure(xtal::Crystal)
	write_vtk(xtal.box, "temp_unit_cell.vtk")
	no_pb = deepcopy(xtal)
	drop_cross_pb_bonds!(no_pb)
	write_mol2(no_pb, filename="temp_view.mol2")
	output = nothing
	try
		output = viewfile("temp_view.mol2", "mol2", vtkcell="temp_unit_cell.vtk")
	catch
		output = viewfile("temp_view.mol2", "mol2", vtkcell="temp_unit_cell.vtk", html=true)
	finally
		rm("temp_unit_cell.vtk")
		rm("temp_view.mol2")
	end
	return output
end


# function to visualize a moiety in the notebook
function view_query_or_replacement(filename::String)
	moty = moiety(filename) # load the moiety
	for i in 1:moty.atoms.n # fix H! atom bonding bug (tagged atom covalent radius too large in JMol)
		if moty.atoms.species[i] == :H!
			moty.atoms.species[i] = :He
		end
		if moty.atoms.species[i] == :C!
			moty.atoms.species[i] = :Ne
		end
	end
	# write temporary modified file, view it, and delete it
	filename = joinpath(PoreMatMod.rc[:paths][:moieties], "temp_" * filename)
	write_xyz(moty, filename)
	output = nothing
	try
		output = viewfile(filename, "xyz")
	catch
		output = viewfile(filename, "xyz", html=true)
	finally
		rm(filename)
	end
	return output
end


# function to print the contents of a moiety file
function display_query_or_replacement_file(filename::String)
	filename = joinpath(PoreMatMod.rc[:paths][:moieties],
		                filename)

	println("contents of: ", filename, "\n")
	open(filename, "r") do io
		print(read(io, String))
	end
end
