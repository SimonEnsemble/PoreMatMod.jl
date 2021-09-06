function view_structure(xtal::Crystal)
	write_vtk(xtal.box, "unit_cell.vtk")
	no_pb = deepcopy(xtal)
	drop_cross_pb_bonds!(no_pb)
	write_mol2(no_pb, filename="view.mol2")
	viewfile("view.mol2", "mol2", vtkcell="unit_cell.vtk")
end

function view_query_or_replacement(filename::String)
	filename = joinpath(PoreMatMod.rc[:paths][:moieties],
		                      filename)
	viewfile(filename, "xyz")
end

function display_query_or_replacement_file(filename::String)
	filename = joinpath(PoreMatMod.rc[:paths][:moieties],
		                filename)

	println("contents of: ", filename, "\n")
	open(filename, "r") do io
		print(read(io, String))
	end
end
