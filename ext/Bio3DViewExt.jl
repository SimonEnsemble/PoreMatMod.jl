module Bio3DViewExt

using PoreMatMod
using Bio3DView: viewfile

include("shared/required_files.jl")

__init__() = check_example_data()

# function to visualize a crystal in the notebook
function view_structure(xtal::Crystal; drop_cross_pb=true)
    # write the box mesh
    write_vtk(xtal.box, "temp_unit_cell.vtk")
    # drop symmetry info and charges
    x = deepcopy(
        Crystal(
            xtal.name,
            xtal.box,
            xtal.atoms,
            Charges{Frac}(0),
            xtal.bonds,
            Xtals.SymmetryInfo()
        )
    )
    # drop cross-boundary bonds (they don't render correctly)
    if drop_cross_pb
        # drop the cross-boundary bonds
        drop_cross_pb_bonds!(x)
    end
    write_mol2(x; filename="temp_view.mol2")
    output = nothing
    try
        output = viewfile("temp_view.mol2", "mol2"; vtkcell="temp_unit_cell.vtk")
    catch
        output = viewfile("temp_view.mol2", "mol2"; vtkcell="temp_unit_cell.vtk", html=true)
    finally
        rm("temp_unit_cell.vtk")
        rm("temp_view.mol2")
    end
    return output
end

# function to visualize a moiety in the notebook
function view_query_or_replacement(filename::String)
    moty = moiety(filename) # load the moiety
    for i in 1:(moty.atoms.n) # fix H! atom bonding bug (tagged atom covalent radius too large in JMol)
        if moty.atoms.species[i] == :H!
            moty.atoms.species[i] = :He
        end
        if moty.atoms.species[i] == :C!
            moty.atoms.species[i] = :Ne
        end
    end
    # write temporary modified file, view it, and delete it
    filename = joinpath(rc[:paths][:moieties], "temp_" * filename)
    write_xyz(moty, filename)
    output = nothing
    try
        output = viewfile(filename, "xyz")
    catch
        output = viewfile(filename, "xyz"; html=true)
    finally
        rm(filename)
    end
    return output
end

export view_query_or_replacement, view_structure

end
