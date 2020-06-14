using PorousMaterials
function moiety(name::String)
    # Generates a moiety (Crystal) from an .xyz file found in path_to_moieties

    box = unit_cube()
    fx = Frac(read_xyz("moieties/$(name).xyz"), box)
    moiety = Crystal(name, box, fx, Charges{Frac}(0))
    infer_bonds!(moiety, false)

    return moiety
end
