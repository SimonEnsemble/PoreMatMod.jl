module temp
using LightGraphs, Test, PorousMaterials, Revise
using MOFun
@testset "bonding" begin
xtal = Crystal("IRMOF-1.cif")
strip_numbers_from_atom_labels!(xtal)
infer_bonds!(xtal, true)
s_moty = moiety("find-replace/2-!-p-phenylene")

nb_bonds(xtal) = ne(xtal.bonds)

# test that a "no-op" leaves the number of bonds unchanged
r_moty = moiety("p-phenylene")
@test nb_bonds((s_moty => r_moty) ∈ xtal) == nb_bonds(xtal)

# test that adding a new moiety increases the number of bonds correctly
r_moty = moiety("2-acetylamido-p-phenylene")
@test ne(((s_moty => r_moty, 1) ∈ xtal).bonds) ==
    (ne(xtal.bonds) - ne(s_moty.bonds) + ne(r_moty.bonds))
end # testset
end # module
