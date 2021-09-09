using Logging

# list of required files for examples
required_files = Dict(:crystals => ["IRMOF-1.cif"],
                      :moieties => ["2-!-p-phenylene.xyz", "2-acetylamido-p-phenylene.xyz"]
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

# function to visualize a crystal in the notebook
function view_structure(xtal::Crystal)
	write_vtk(xtal.box, "temp_unit_cell.vtk")
	no_pb = deepcopy(xtal)
	drop_cross_pb_bonds!(no_pb)
	write_mol2(no_pb, filename="temp_view.mol2")
	output = viewfile("temp_view.mol2", "mol2", vtkcell="temp_unit_cell.vtk")
	rm("temp_unit_cell.vtk")
	rm("temp_view.mol2")
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
	output = viewfile(filename, "xyz")
	rm(filename)
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
