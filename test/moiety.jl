module Moiety_Test

using Test, PoreMatMod

# function for a specific test case. NOT a generally useful function!
# frag and moty inputs are expected to be permuted sets of atoms with unique species
function test_is_equiv(frag::Crystal, moty::Crystal)::Bool
    for (f, f_label) in enumerate(frag.atoms.species)
        for (m, m_label) in enumerate(moty.atoms.species)
            if f_label == m_label
                for i in 1:3
                    if frag.atoms.coords.xf[i,f] ≠ moty.atoms.coords.xf[i,m]
                        return false # detected an atom w/ inconsistent coordinates
                    end
                end
            end
        end
    end
    return true
end


@testset "Moiety Tests" begin
    moiety_bcfm = moiety("S-bromochlorofluoromethane.xyz")
    moiety_2!bcfm = moiety("!-S-bromochlorofluoromethane.xyz")
    @test moiety_bcfm.atoms.species == [:C, :Cl, :F, :Br, :H]
    @test moiety_2!bcfm.atoms.species == [:C, :Cl, :F, :Br, :H!]
end

end # module
