module MOFun_Test

using Test
using PorousMaterials, MOFun, Moiety

@testset "MOFun Tests" begin
    irmof1 = Crystal("IRMOF-1.cif")
    strip_numbers_from_atom_labels!(irmof1)
    infer_bonds!(irmof1, true)
    timil125 = Crystal("Ti-MIL-125.cif")
    strip_numbers_from_atom_labels!(timil125)
    infer_bonds!(timil125, true)
    moty = moiety("p-phenylene")
    moty_w_R_grp = moiety("find-replace/2-!-p-phenylene")
    results1 = moty ∈ irmof1
    @test length(results1) == 96
    @test length([result for result in results1 if result.configuration.location == 1]) == 4
    results2 = moty ∈ timil125
    @test length(results2) == 48
    @test length([result for result in results2 if result.configuration.location == 1]) == 4
    results3 = moty_w_R_grp ∈ timil125
    @test length(results3) == 48
    @test length([result for result in results3 if result.configuration.location == 1]) == 4
end
end
