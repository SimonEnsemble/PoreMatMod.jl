module Moiety_Test

using Test, PorousMaterials, LightGraphs
include("../src/moiety.jl")

@testset "Moiety Tests" begin
    @test moiety("find-replace/2-!-p-phenylene") ≠ nothing
    @test moiety("p-phenylene") ≠ nothing
end

end
