module Moiety_Test

using Test, Revise
using PorousMaterials, LightGraphs, MOFun
include("../src/moiety.jl")


# function for a specific test case. NOT a generally useful function!
function test_is_equiv(frag::Crystal, moty::Crystal)::Bool
    @debug "fragment" frag.atoms.species frag.atoms.coords
    @debug "moiety" moty.atoms.species moty.atoms.coords
    for (f, f_label) in enumerate(frag.atoms.species)
        for (m, m_label) in enumerate(moty.atoms.species)
            if f_label == m_label
                @debug "testing" f_label m_label frag.atoms.coords[f] moty.atoms.coords[m]
                for i in 1:3
                    if frag.atoms.coords.xf[i,f] ≠ moty.atoms.coords.xf[i,m]
                        return false
                    end
                end
            end
        end
    end
    return true
end


@testset "Moiety Tests" begin
    fragment_bcfm = read_fragment_from_xyz("test/S-bromochlorofluoromethane", fragment_location=PATH_TO_MOIETIES)
    moiety_bcfm = moiety("test/S-bromochlorofluoromethane")
    fragment_2!bcfm = read_fragment_from_xyz("test/!-S-bromochlorofluoromethane", fragment_location=PATH_TO_MOIETIES)
    moiety_2!bcfm = moiety("test/!-S-bromochlorofluoromethane")
    # check that the species/coordinates correspondence is identical between fragments and moieties
    @test test_is_equiv(fragment_bcfm, moiety_bcfm)
    @test test_is_equiv(fragment_2!bcfm, moiety_2!bcfm)
end

end # module
