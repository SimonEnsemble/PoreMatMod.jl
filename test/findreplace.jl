module MOFun_Test

using LightGraphs, Test, Revise
using Xtals, MOFun

@testset "substructure_search" begin
irmof1 = Crystal("IRMOF-1.cif")
strip_numbers_from_atom_labels!(irmof1)
infer_bonds!(irmof1, true)
timil125 = Crystal("Ti-MIL-125.cif")
strip_numbers_from_atom_labels!(timil125)
infer_bonds!(timil125, true)
p_phenylene = moiety("p-phenylene")
p_phenylene_w_R_grp = moiety("2-!-p-phenylene")
search1 = p_phenylene ∈ irmof1
@test nb_isomorphisms(search1) == 96
@test nb_locations(search1) == 24
@test nb_configs_at_loc(search1)[1] == 4
@test search1.results[1].isomorphism[1] ==
    [233, 306, 318, 245, 185, 197, 414, 329, 402, 341]
search2 = p_phenylene ∈ timil125
@test nb_isomorphisms(search2) == 48
@test nb_locations(search2) == 12
@test nb_configs_at_loc(search2)[1] == 4
@test search2.results[1].isomorphism[1] ==
    [8, 140, 144, 141, 7, 133, 186, 185, 190, 189]
search3 = p_phenylene_w_R_grp ∈ timil125
@test nb_isomorphisms(search3) == 48
@test nb_locations(search3) == 12
@test nb_configs_at_loc(search3)[1] == 4
@test search3.results[1].isomorphism[1] ==
    [8, 140, 144, 141, 7, 133, 185, 190, 189, 186]
s_moty = moiety("!-S-bromochlorofluoromethane")
parent = moiety("S-bromochlorofluoromethane")
search = s_moty ∈ parent
@test search.results[1].isomorphism[1] == [1, 4, 2, 3, 5]
@test parent.atoms.species[search.results[1].isomorphism[1]][1:4] ==
    s_moty.atoms.species[1:4] &&
    s_moty.atoms.species[5] == :H! &&
    parent.atoms.species[search.results[1].isomorphism[1]][5] == :H
s_moty = moiety("p-phenylene")
parent = Crystal("IRMOF-1_one_ring.cif")
strip_numbers_from_atom_labels!(parent)
infer_bonds!(parent, true)
search = s_moty ∈ parent
@test search.results[1].isomorphism[1] == [34, 38, 39, 36, 26, 27, 51, 46, 50, 48]
end # test set: substructure_search

@testset "find_repalce" begin
parent = Crystal("IRMOF-1.cif")
strip_numbers_from_atom_labels!(parent)
infer_bonds!(parent, true)
s_moty = moiety("2-!-p-phenylene")
r_moty = moiety("2-acetylamido-p-phenylene")
new_xtal = (s_moty => r_moty) ∈ parent
@test new_xtal.atoms.n == 592
new_xtal = (s_moty => r_moty, 1) ∈ parent
@test new_xtal.atoms.n == 431
new_xtal = (s_moty => r_moty, [2,3]) ∈ parent
@test new_xtal.atoms.n == 438
new_xtal = (s_moty => r_moty, [2,3,4], [1,1,1]) ∈ parent
@test new_xtal.atoms.n == 445
r_moty = moiety("p-phenylene")
new_xtal = (s_moty => r_moty) ∈ parent
@test ne(new_xtal.bonds) == ne(parent.bonds)
r_moty = moiety("2-acetylamido-p-phenylene")
new_xtal = (s_moty => r_moty, 1) ∈ parent
@test ne(new_xtal.bonds) == (ne(parent.bonds) - ne(s_moty.bonds) + ne(r_moty.bonds))
xtal = Crystal("IRMOF-1.cif")
strip_numbers_from_atom_labels!(xtal)
infer_bonds!(xtal, true)
s_moty = moiety("2-!-p-phenylene")
nb_bonds(xtal) = ne(xtal.bonds)
# test that a "no-op" leaves the number of bonds unchanged
r_moty = moiety("p-phenylene")
@test nb_bonds((s_moty => r_moty) ∈ xtal) == nb_bonds(xtal)
# test that adding a new moiety increases the number of bonds correctly
r_moty = moiety("2-acetylamido-p-phenylene")
@test ne(((s_moty => r_moty, 1) ∈ xtal).bonds) ==
    (ne(xtal.bonds) - ne(s_moty.bonds) + ne(r_moty.bonds))
end # test set: find_repalce

end # module
