### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ d757a3a8-1ee3-4afd-b0bd-9d8429250f96
begin
    import Pkg
    Pkg.add(url="https://github.com/SimonEnsemble/PoreMatMod.jl")
end

# ╔═╡ 322d1bd1-86ae-47e7-b49d-38687eb3885a
using PoreMatMod

# ╔═╡ 494d2664-b5cf-43fe-850b-538a090dd0e8
using PoreMatMod.ExampleHelpers

# ╔═╡ 095d18b8-b8ab-4bd8-89dc-4382f9375b42
begin
    parent = Crystal("UiO-66.cif")
    for x in parent.atoms.coords.xf
        if x == 1.
            x -= eps(Float64)
        elseif x == 0.
            x += eps(Float64)
        end
    end
end

# ╔═╡ ac94ddb2-314c-49a7-b0bb-4371f93c0429
infer_bonds!(parent, true)

# ╔═╡ ecb86319-1c5c-48f8-98b9-8bccca5e9953
SBU = moiety("SBU.xyz")

# ╔═╡ d7b5abf5-13fb-4644-9502-6c54106be51a
bimetallic_SBU = moiety("bimetallic_SBU.xyz")

# ╔═╡ 44b0cc75-036b-468f-9da7-04c2d97536a9
begin
    multi_metal = replace(parent, SBU => bimetallic_SBU)
    translate_by!(multi_metal.atoms.coords, Frac([-0.25, -0.25, 0.]))
end

# ╔═╡ 780437f4-3a42-4aa2-b2e8-599d8d232dc2
write_cif(multi_metal, "multi-metal.cif")

# ╔═╡ ad888be3-2123-4d2f-865a-ff7a63cebf9e
begin
    shifted_uio66 = deepcopy(parent)
    translate_by!(shifted_uio66.atoms.coords, Frac([-0.25, -0.25, 0.]))
    write_cif(shifted_uio66, "shifted_parent.cif")
end

# ╔═╡ Cell order:
# ╠═d757a3a8-1ee3-4afd-b0bd-9d8429250f96
# ╠═322d1bd1-86ae-47e7-b49d-38687eb3885a
# ╠═494d2664-b5cf-43fe-850b-538a090dd0e8
# ╠═095d18b8-b8ab-4bd8-89dc-4382f9375b42
# ╠═ac94ddb2-314c-49a7-b0bb-4371f93c0429
# ╠═ecb86319-1c5c-48f8-98b9-8bccca5e9953
# ╠═d7b5abf5-13fb-4644-9502-6c54106be51a
# ╠═44b0cc75-036b-468f-9da7-04c2d97536a9
# ╠═780437f4-3a42-4aa2-b2e8-599d8d232dc2
# ╠═ad888be3-2123-4d2f-865a-ff7a63cebf9e
