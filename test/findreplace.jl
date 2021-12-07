module PoreMatMod_Test

using Test, Graphs, PoreMatMod, LinearAlgebra

@testset "non-p1 symmetry" begin
    parent = Crystal("NiPyC_fragment_trouble.cif", convert_to_p1=false)
    infer_bonds!(parent, true)

    query = moiety("PyC.xyz")
    replacement = moiety("PyC-CH3.xyz")
    prim_child = replace(parent, query => replacement)

    write_cif(prim_child, joinpath(rc[:paths][:crystals], "prim_child.cif"))
    child = Crystal("prim_child.cif")

    @test child.atoms.n == 33
    @test ne(prim_child.bonds) == 34
    @test   prim_child.symmetry.operations == parent.symmetry.operations && 
            prim_child.symmetry.space_group == parent.symmetry.space_group && 
            prim_child.symmetry.is_p1 == parent.symmetry.is_p1
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
@test search1.isomorphisms[1][1] == [233, 306, 318, 245, 185, 197, 414, 329, 402, 341]
search2 = p_phenylene ∈ timil125
@test nb_isomorphisms(search2) == 48
@test nb_locations(search2) == 12
@test nb_ori_at_loc(search2)[1] == 4
@test search2.isomorphisms[1][1] == [8, 140, 144, 141, 7, 133, 186, 185, 190, 189]
search3 = p_phenylene_w_R_grp ∈ timil125
@test nb_isomorphisms(search3) == 48
@test nb_locations(search3) == 12
@test nb_ori_at_loc(search3)[1] == 4
@test search3.isomorphisms[1][1] == [8, 140, 144, 141, 7, 133, 185, 190, 189, 186]
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

@testset "platform-specific output consistency" begin
"""
        !!! DEVELOPER NOTE !!!

    This test is weird. The exact coordinates that come out of the find/replace operation being tested depend on the system in use.
    Specifically, Intel Haswell CPUs give slightly different results than Intel SkylakeX. Why? Not sure. They are visually indiscernible.
    So, different manually-authenticated outputs are used for the test.
    If the machine is running OpenBLAS 0.3.10 and LAPACK 3.9.0 on a SkylakeX CPU, the test is vs. skylakex_acetamido_IRMOF-1.cif.
    If the machine is running those OpenBLAS/LAPACK versions on a Haswell CPU, the test is vs. haswell_acetamido_IRMOF-1.cif.
    The Haswell coords are used as a fallback on machines matching neither description, and the test's results are displayed but ignored.
"""

environment = :unknown

try
    libblas = LinearAlgebra.BLAS.libblas
    blas_string = LinearAlgebra.BLAS.openblas_get_config()
    liblapack = LinearAlgebra.LAPACK.liblapack
    lapack_version = LinearAlgebra.LAPACK.version()

    if libblas == liblapack == "libopenblas64_" && contains(blas_string, "OpenBLAS") && contains(blas_string, "0.3.10") && lapack_version == v"3.9.0"
        if contains(blas_string, "SkylakeX")
            environment = :SkylakeX
        elseif contains(blas_string, "Haswell")
            environment = :Haswell
        end
    end
catch
end

if environment == :unknown
    @warn "The present environment has undefined behavior for this test. The test will run and its results will be displayed, but ignored. The user should visually inspect test_acetamido_IRMOF-1.cif"
end

# test that the coordinates resulting from a specific replacement are the same as a verified test run
parent = Crystal("IRMOF-1.cif")
infer_bonds!(parent, true)
query = moiety("2-!-p-phenylene.xyz")
replacement = moiety("2-acetylamido-p-phenylene.xyz")
xtal1 = replace(parent, query => replacement, loc=[2,4,6,8], ori=[1,2,3,4])
write_cif(xtal1, "test_acetamido_IRMOF-1.cif") # for CI artifact collection

if environment == :SkylakeX
    xtal2 = Crystal("skylake_acetamido_IRMOF-1.cif")
else
    xtal2 = Crystal("haswell_acetamido_IRMOF-1.cif")
end

result = all(isapprox.(xtal1.atoms.coords.xf, xtal2.atoms.coords.xf, atol=0.001))

if environment ≠ :unknown
    @test result
elseif result
    @info "Test passed."
else
    @warn "Test failed on unknown environment. See github.com/SimonEnsemble/PoreMatMod.jl/issues/92 for more info."
end
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
