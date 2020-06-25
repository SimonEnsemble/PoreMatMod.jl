# moiety.jl
# Adrian Henle, Cory Simon 2020

using PorousMaterials

## TODO
"""
    # add path_to_moieties to PorousMaterials
    # moiety() tests:
        # generate vtk
        # generate manual graph representation
        # verify that the two are identical
"""
## #


"""
    Generates a moiety (Crystal) from an .xyz file found in path_to_moieties
"""
function moiety(name::String)

    box = unit_cube()
    fx = Frac(read_xyz("$(name).xyz"), box)
    moiety = Crystal(name, box, fx, Charges{Frac}(0))
    infer_bonds!(moiety, false)

    return moiety
end
