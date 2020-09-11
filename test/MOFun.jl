module MOFun_Test

using Test, Revise
using PorousMaterials, MOFun

@testset "substructure_search" begin

irmof1 = Crystal("IRMOF-1.cif")
strip_numbers_from_atom_labels!(irmof1)
infer_bonds!(irmof1, true)
timil125 = Crystal("Ti-MIL-125.cif")
strip_numbers_from_atom_labels!(timil125)
infer_bonds!(timil125, true)
p_phenylene = moiety("p-phenylene")
p_phenylene_w_R_grp = moiety("find-replace/2-!-p-phenylene")
search1 = p_phenylene ∈ irmof1
@test search1.num_isomorphisms == 96
@test search1.num_locations == 24
@test search1.num_orientations[1] == 4
@test search1.results[1].isomorphism ==
    [259, 284, 295, 272, 211, 220, 391, 355, 380, 368]
search2 = p_phenylene ∈ timil125
@test search2.num_isomorphisms == 48
@test search2.num_locations == 12
@test search2.num_orientations[1] == 4
@test search2.results[1].isomorphism ==
    [155, 170, 162, 179, 59, 66, 207, 208, 205, 206]
search3 = p_phenylene_w_R_grp ∈ timil125
@test search3.num_isomorphisms == 48
@test search3.num_locations == 12
@test search3.num_orientations[1] == 4
@test search3.results[1].isomorphism ==
    [155, 170, 162, 179, 59, 66, 208, 205, 206, 207]
s_moty = moiety("test/!-S-bromochlorofluoromethane")
parent = moiety("test/S-bromochlorofluoromethane")
search = s_moty ∈ parent
@test search.results[1].isomorphism == [1, 3, 4, 5, 2]
@test parent.atoms.species[search.results[1].isomorphism][1:4] ==
    s_moty.atoms.species[1:4] &&
    s_moty.atoms.species[5] == :H! &&
    parent.atoms.species[search.results[1].isomorphism][5] == :H
s_moty = moiety("p-phenylene")
parent = Crystal("IRMOF-1_one_ring.cif")
strip_numbers_from_atom_labels!(parent)
infer_bonds!(parent, true)
search = s_moty ∈ parent
@test search.results[1].isomorphism == [34, 38, 39, 36, 26, 27, 51, 46, 50, 48]

end # test set: substructure_search

@testset "find_repalce" begin
    @test true
end # test set: find_repalce

end # module
