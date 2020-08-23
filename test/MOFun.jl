module MOFun_Test

using Test, Revise
using PorousMaterials, MOFun

@testset "MOFun Tests" begin
    irmof1 = Crystal("IRMOF-1.cif")
    strip_numbers_from_atom_labels!(irmof1)
    infer_bonds!(irmof1, true)
    timil125 = Crystal("Ti-MIL-125.cif")
    strip_numbers_from_atom_labels!(timil125)
    infer_bonds!(timil125, true)
    p_phenylene = moiety("p-phenylene")
    p_phenylene_w_R_grp = moiety("find-replace/2-!-p-phenylene")
    search1 = p_phenylene ∈ irmof1
    @test length(search1.results) == 96
    @test length([result for result in search1.results if result.configuration.location == 1]) == 4
    search2 = p_phenylene ∈ timil125
    @test length(search2.results) == 48
    @test length([result for result in search2.results if result.configuration.location == 1]) == 4
    search3 = p_phenylene_w_R_grp ∈ timil125
    @test length(search3.results) == 48
    @test length([result for result in search3.results if result.configuration.location == 1]) == 4
end
end
