module PoreMatMod_Test

using Test, PoreMatMod

@testset "correct_missing_Hs" begin
    include("../examples/correct_missing_Hs.jl")
    @test true
end

@testset "disorder_and_guests" begin
    include("../examples/disorder_and_guests.jl")
    @test true
end

@testset "make_hypothetical_MOF" begin
    include("../examples/make_hypothetical_MOF.jl")
    @test true
end

@testset "missing_linker_defect" begin
    include("../examples/missing_linker_defect.jl")
    @test true
end

@testset "replacement_modes" begin
    include("../examples/replacement_modes.jl")
    @test true
end

end
