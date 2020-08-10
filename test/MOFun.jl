module MOFun_Test

using MOFun, PorousMaterials, Test

irmof1 = Crystal("IRMOF-1.cif")
strip_numbers_from_atom_labels!(irmof1)
## BUG
# Getting bonds with this line raises an exception:
# infer_geometry_based_bonds!(irmof1, true)
infer_bonds!(irmof1, true)
timil125 = Crystal("Ti-MIL-125.cif")
strip_numbers_from_atom_labels!(timil125)
infer_bonds!(timil125, true)
moty = moiety("p-phenylene")

@testset "MOFun Tests" begin
    @test length(substructure_search(moty, irmof1)) == 96
    @test length(substructure_search(moty, timil125)) == 48
end
end
