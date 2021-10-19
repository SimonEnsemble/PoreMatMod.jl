module PoreMatMod_Test

using Test, LightGraphs, PoreMatMod

@testset "substructure_search" begin
irmof1 = Crystal("IRMOF-1.cif")
strip_numbers_from_atom_labels!(irmof1)
infer_bonds!(irmof1, true)
timil125 = Crystal("Ti-MIL-125.cif")
strip_numbers_from_atom_labels!(timil125)
infer_bonds!(timil125, true)
p_phenylene = moiety("p-phenylene.xyz")
p_phenylene_w_R_grp = moiety("2-!-p-phenylene.xyz")
search1 = p_phenylene ∈ irmof1
@test nb_isomorphisms(search1) == 96
@test nb_locations(search1) == 24
@test nb_ori_at_loc(search1)[1] == 4
@test search1.isomorphisms[1][1] ==
    [233, 306, 318, 245, 185, 197, 414, 329, 402, 341]
search2 = p_phenylene ∈ timil125
@test nb_isomorphisms(search2) == 48
@test nb_locations(search2) == 12
@test nb_ori_at_loc(search2)[1] == 4
@test search2.isomorphisms[1][1] ==
    [8, 140, 144, 141, 7, 133, 186, 185, 190, 189]
search3 = p_phenylene_w_R_grp ∈ timil125
@test nb_isomorphisms(search3) == 48
@test nb_locations(search3) == 12
@test nb_ori_at_loc(search3)[1] == 4
@test search3.isomorphisms[1][1] ==
    [8, 140, 144, 141, 7, 133, 185, 190, 189, 186]
query = moiety("!-S-bromochlorofluoromethane.xyz")
parent = moiety("S-bromochlorofluoromethane.xyz")
search = query ∈ parent
@test search.isomorphisms[1][1] == [1, 2, 3, 4, 5]
@test parent.atoms.species[search.isomorphisms[1][1]][1:4] ==
    query.atoms.species[1:4] &&
    query.atoms.species[5] == :H! &&
    parent.atoms.species[search.isomorphisms[1][1]][5] == :H
query = moiety("p-phenylene.xyz")
parent = Crystal("IRMOF-1_one_ring.cif")
strip_numbers_from_atom_labels!(parent)
infer_bonds!(parent, true)
search = query ∈ parent
@test search.isomorphisms[1][1] == [34, 38, 39, 36, 26, 27, 51, 46, 50, 48]
end # test set: substructure_search

@testset "find_repalce" begin
parent = Crystal("IRMOF-1.cif")
strip_numbers_from_atom_labels!(parent)
infer_bonds!(parent, true)
query = moiety("2-!-p-phenylene.xyz")
replacement = moiety("2-acetylamido-p-phenylene.xyz")
new_xtal = replace(parent, query => replacement)
@test new_xtal.atoms.n == 592
new_xtal = replace(parent, query => replacement, nb_loc=1)
@test new_xtal.atoms.n == 431
new_xtal = replace(parent, query => replacement, loc=[2,3])
@test new_xtal.atoms.n == 438
new_xtal = replace(parent, query => replacement, loc=[2,3,4], ori=[1,1,1])
@test new_xtal.atoms.n == 445
replacement = moiety("p-phenylene.xyz")
new_xtal = replace(parent, query => replacement)
@test ne(new_xtal.bonds) == ne(parent.bonds)
replacement = moiety("2-acetylamido-p-phenylene.xyz")
new_xtal = replace(parent, query => replacement, nb_loc=1)
@test ne(new_xtal.bonds) == (ne(parent.bonds) - ne(query.bonds) + ne(replacement.bonds))
xtal = Crystal("IRMOF-1.cif")
strip_numbers_from_atom_labels!(xtal)
infer_bonds!(xtal, true)
query = moiety("2-!-p-phenylene.xyz")
nb_bonds(xtal) = ne(xtal.bonds)
# test that a "no-op" leaves the number of bonds unchanged
replacement = moiety("p-phenylene.xyz")
@test nb_bonds(replace(xtal, query => replacement)) == nb_bonds(xtal)
# test that adding a new moiety increases the number of bonds correctly
replacement = moiety("2-acetylamido-p-phenylene.xyz")
@test ne((replace(xtal, query => replacement, nb_loc=1)).bonds) ==
    (ne(xtal.bonds) - ne(query.bonds) + ne(replacement.bonds))
# test that the coordinates resulting from a specific replacement are the same as a verified test run
parent = Crystal("IRMOF-1.cif")
infer_bonds!(parent, true)
query = moiety("2-!-p-phenylene.xyz")
replacement = moiety("2-acetylamido-p-phenylene.xyz")
xtal1 = replace(parent, query => replacement, loc=[2,4,6,8], ori=[1,2,3,4])
xtal2 = Crystal("verified_acetamido_IRMOF-1.cif")
write_cif(xtal1, "test_acetamido_IRMOF-1.cif") ## TODO remove (added for debugging)
@test all(isapprox.(xtal1.atoms.coords.xf, xtal2.atoms.coords.xf, rtol=1e-3))
end # test set: find_repalce

end # module
