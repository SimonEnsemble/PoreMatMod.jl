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
    @test length(moty ∈ irmof1) == 24
    @test length(moty ∈ timil125) == 12
    @test length(moty_w_R_grp ∈ timil125) == 12
end
end
