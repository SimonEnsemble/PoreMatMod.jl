# moiety.jl
# Adrian Henle, Cory Simon 2020

using PorousMaterials
PATH_TO_MOIETIES = joinpath(pwd(), "data/moieties")

## TODO
"""
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
    fx = Frac(read_xyz(joinpath(PATH_TO_MOIETIES, "$(name).xyz")), box)
    moiety = Crystal(name, box, fx, Charges{Frac}(0))
    infer_bonds!(moiety, false)

    return moiety
end
