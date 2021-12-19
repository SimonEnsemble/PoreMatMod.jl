module PoreMatMod_Test

using Test, Graphs, PoreMatMod, LinearAlgebra

@testset "replacement split across PB" begin
    parent = Crystal("NiPyC_fragment_trouble.cif")
    query = moiety("PyC.xyz")
    replacement = moiety("PyC-CH3.xyz")

    @test_throws AssertionError replace(parent, query => replacement)
    infer_bonds!(parent, true)
    child = replace(parent, query => replacement)
    
    @test ne(child.bonds) == ne(parent.bonds) + 3 * 2 # added H -> CH3 on two PyC ligands

    # test conglomerate worked.
    remove_bonds!(child)
    infer_bonds!(child, true)
    @test ne(child.bonds) == ne(parent.bonds) + 3 * 2 # added H -> CH3 on two PyC ligands
end

@testset "split across PB, non-connected replacement" begin
    parent = Crystal("NiPyC_fragment_trouble.cif", convert_to_p1=false)
    infer_bonds!(parent, true)
    
    query = moiety("PyC_split.xyz")
    replacement = moiety("PyC_split_replacement.xyz")
    child = replace(parent, query => replacement)
    
    # ensure O atoms aligned
    parent_Os = parent[parent.atoms.species .== :O]
    child_Os  = child[child.atoms.species .== :O]
    nb_Os = parent_Os.atoms.n
    @test parent_Os.atoms.n == nb_Os
    
    min_distances = [Inf for _ = 1:nb_Os] # from parent
    for i = 1:nb_Os
        for j = 1:nb_Os
            r = norm(parent_Os.atoms.coords.xf[:, i] - child_Os.atoms.coords.xf[:, j])
            if r < min_distances[i] 
                min_distances[i] = r
            end
        end
    end
    @test all(min_distances .< 0.01)
end

@testset "replacement spans twice the unit cell" begin
    parent = Crystal("NiPyC_fragment_trouble.cif")
    infer_bonds!(parent, true)

    query = moiety("PyC.xyz")
    replacement = moiety("PyC-long_chain.xyz")
    
    @test_throws Exception replace(parent, query => replacement)
end

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
    
    @test search1.isomorphisms[1][1] == Dict([q => p for (q, p) in enumerate([233, 306, 318, 245, 185, 197, 414, 329, 402, 341])])
    
    search2 = p_phenylene ∈ timil125
    
    @test nb_isomorphisms(search2) == 48
    
    @test nb_locations(search2) == 12
    
    @test nb_ori_at_loc(search2)[1] == 4
    
    @test search2.isomorphisms[1][1] == Dict([q => p for (q, p) in enumerate([8, 140, 144, 141, 7, 133, 186, 185, 190, 189])])
    
    search3 = p_phenylene_w_R_grp ∈ timil125
    
    @test nb_isomorphisms(search3) == 48
    
    @test nb_locations(search3) == 12
    
    @test nb_ori_at_loc(search3)[1] == 4
    
    @test search3.isomorphisms[1][1] == Dict([q => p for (q, p) in enumerate([8, 140, 144, 141, 7, 133, 185, 190, 189, 186])])
    
    query = moiety("!-S-bromochlorofluoromethane.xyz")
    parent = moiety("S-bromochlorofluoromethane.xyz")
    search = query ∈ parent
    
    @test search.isomorphisms[1][1] == Dict([q => p for (q, p) in enumerate([1, 2, 3, 4, 5])])
    
    @test [parent.atoms.species[search.isomorphisms[1][1][i]] for i in 1:4] == query.atoms.species[1:4] &&
        query.atoms.species[5] == :H! && parent.atoms.species[search.isomorphisms[1][1][5]] == :H
    
    query = moiety("p-phenylene.xyz")
    parent = Crystal("IRMOF-1_one_ring.cif")
    strip_numbers_from_atom_labels!(parent)
    infer_bonds!(parent, true)
    search = query ∈ parent
    
    @test search.isomorphisms[1][1] == Dict([q => p for (q, p) in enumerate([34, 38, 39, 36, 26, 27, 51, 46, 50, 48])])
end # test set: substructure_search


@testset "find_replace" begin
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
end

@testset "conglomerate test" begin
    xtal = Crystal("conglomerate_test.cif")
    infer_bonds!(xtal, false)
    @test ne(xtal.bonds) == 1
    remove_bonds!(xtal)
    @test ne(xtal.bonds) == 0
    infer_bonds!(xtal, true)
    @test ne(xtal.bonds) == 5
    PoreMatMod.conglomerate!(xtal)
    remove_bonds!(xtal)
    translate_by!(xtal.atoms.coords, Frac([0.5, 0.5, 0.5]))
    infer_bonds!(xtal, false)
    @test ne(xtal.bonds) == 5
end


@testset "remove duplicates" begin
    parent = moiety("ADC.xyz")
    new_box = replicate(unit_cube(), (10, 10, 10))
    new_atoms = Frac(Cart(parent.atoms, parent.box), new_box)
    new_charges = Frac(Cart(parent.charges, parent.box), new_box)
    parent = Crystal(parent.name, new_box, new_atoms, new_charges)
    infer_bonds!(parent, false)
    query = moiety("naphthyl_fragment.xyz")
    replacement = moiety("F_naphthyl_fragment.xyz")
    child = replace(parent, query => replacement, nb_loc=2, remove_duplicates=true, reinfer_bonds=true)

    @test child.atoms.n == parent.atoms.n

    @test nv(child.bonds) == nv(parent.bonds) && ne(child.bonds) == ne(parent.bonds)
end

end
