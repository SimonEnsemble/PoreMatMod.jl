# load the data if missing
if !isdir(joinpath(pwd(), "data"))
	mkpath(joinpath(pwd(), "data", "crystals"))
	mkpath(joinpath(pwd(), "data", "moieties"))
	using ZipFile
	download("https://github.com/SimonEnsemble/PoreMatMod.jl/raw/0093f074fca8c14515ef3125f5fb03b22a0a8303/docs/src/assets/examples/examples.zip", "temp_examples.zip")
	r = ZipFile.Reader("temp_examples.zip")
	for f in r.files
		try
			open(f.name, "w") do io
				write(io, read(f, String))
			end
		catch
		end
	end
	close(r)
	rm("temp_examples.zip")
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
