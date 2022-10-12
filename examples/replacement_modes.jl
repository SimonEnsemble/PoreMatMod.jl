### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 359a6c00-c9a6-441d-b258-55bfb5deb4b5
begin
    import Pkg
    Pkg.add(; url="https://github.com/SimonEnsemble/PoreMatMod.jl")
end

# ╔═╡ a948f8b3-4ec5-40b9-b2c1-fcf5b8ad67fa
# load required packages (Pluto.jl will automatically install them)
using PoreMatMod, PlutoUI

# ╔═╡ 77c63e60-2708-4fef-a6ea-cb394f114b88
begin
    using PoreMatMod.ExampleHelpers
    check_example_data()
end

# ╔═╡ 8d523993-6e85-443a-9949-12030552b457
md"""
## example: generating hypothetical MOF structures with different replacement styles
"""

# ╔═╡ 74b45651-21d5-4332-a4a5-866ea1bb02b8
input_file_message()

# ╔═╡ 3ff0cd63-8bdf-4ba5-92b6-e2a52f7573c2
xtal_folder_message()

# ╔═╡ ccc7fc86-5f48-4cc9-bd5c-349fc1d55b0c
rc[:paths][:crystals]

# ╔═╡ b66dd0d5-5bac-4b23-a33f-17305b0d4728
moiety_folder_message()

# ╔═╡ a76e79ee-c82c-4771-818d-5380a5fb4c18
rc[:paths][:moieties]

# ╔═╡ be86f07e-0669-46c2-8c79-80e65dfcc6f2
md"""
!!! example \"the task\"
    we have the IRMOF-1 crystal structure, and wish to explore appending the acetamido functional group its BDC linkers using different replacement styles.

**Parent crystal structure**: first, we read in the `.cif` file describing the parent IRMOF-1 crystal structure.
"""

# ╔═╡ b53e7c38-d8f5-4f28-a9dc-f7902a86fdb2
begin
    # read in the parent xtal
    parent = Crystal("IRMOF-1.cif") # load .cif file
    infer_bonds!(parent, true)      # infer bonds
    view_structure(parent)          # view structure
end

# ╔═╡ b402e79a-784b-4d8b-82f1-df4fe1cedad1
fragment_construction_note()

# ╔═╡ 9b2d534a-4f78-4950-add1-9ba645669bb9
md"""
**Query fragment**: next, we construct a query fragment as a _p_-phenylene fragment to match that on the BCD linker of the IRMOF-1 parent structure. we mark one hydrogen atom on the query fragment as "masked" (shown in pink) by tagging its species label with `!` in the input. we need to mask this hydrogen atom because it will eventually be replaced by the acetamido functional group.
"""

# ╔═╡ 4be03110-61ab-4cd6-b3fe-7d51ac5ee771
query = moiety("2-!-p-phenylene.xyz");

# ╔═╡ e5eaaef4-13a0-48bd-9f2b-5040b2d10ac1
view_query_or_replacement("2-!-p-phenylene.xyz")

# ╔═╡ 559ef3c3-a176-4d65-8958-810c9b0b32c5
with_terminal() do
    return display_query_or_replacement_file("2-!-p-phenylene.xyz")
end

# ╔═╡ d1aa8a19-3702-40f6-b73e-b9ebfc1a7a71
md"""
**Replacement fragment**: next, we construct a replacement fragment as a modified version of the query fragment (with acetamido group in place of one hydrogen atom).
"""

# ╔═╡ 44e7e665-ae68-4cd4-b45f-138a0fb8910e
replacement = moiety("2-nitro-p-phenylene.xyz");

# ╔═╡ af606f02-fbd5-4033-a5f3-56a8f740b1a1
view_query_or_replacement("2-nitro-p-phenylene.xyz")

# ╔═╡ c8f8dbd3-4191-460d-944f-f5e456ce8b83
with_terminal() do
    return display_query_or_replacement_file("2-nitro-p-phenylene.xyz")
end

# ╔═╡ 127a2bef-3903-4801-bc75-00a6dde2bc6e
md"""
### Replacement Modes

There are several options when performing a replacement that control how the replacement fragments will be positioned in the parent structure.

With all three file inputs loaded (IRMOF-1 as `parent`, 2-!-*p*-phenylene as `query`, and 2-nitro-*p*-phenylene as `replacement`) and a `search` performed, replacements can be made.

`PoreMatMod.jl` has several replacement modes, one of which must be specified.
"""

# ╔═╡ 0e3f7334-9e7f-483b-a0ac-70777902bf51
md"""
!!! note \"Note\"
    Multiple replacements can be done with a single search.
"""

# ╔═╡ b7b46022-aee0-4d51-8a5e-0c8c005f341a
search = query ∈ parent

# ╔═╡ 68557426-3204-45a1-8801-84c459e5ac66
md"""
#### Default: optimal orientation at all locations

Optimal configurations will be chosen for each location in `search.isomorphisms`, so that each occurrence of the `query` in the `parent` is replaced with minimal perturbation of the conserved atoms from the parent structure.
"""

# ╔═╡ 74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
begin
    local child = substructure_replace(search, replacement)
    view_structure(child)
end

# ╔═╡ bed9d62f-d61d-4fb1-9a12-188a864fff21
md"""
#### Optimal Orientation at *n* random locations

If only some of the linkers should be replaced, the `nb_loc` argument lets us specify how many.
"""

# ╔═╡ 9dbf7859-762d-4edc-8300-ac3bba151a8a
begin
    local child = substructure_replace(search, replacement; nb_loc=8)
    view_structure(child)
end

# ╔═╡ 95215fe9-a7e6-4563-94b6-5ed1d5dde03b
md"""
#### Optimal orientation at specific locations

Specific locations are chosen by providing the `loc` argument.
"""

# ╔═╡ 8c22f94c-b235-40b2-8fc8-6052c31a6b6e
begin
    local child = substructure_replace(search, replacement; loc=[17, 18, 19, 20])
    view_structure(child)
end

# ╔═╡ 6297833a-ea26-4958-ad56-fc76aeabee69
md"""
#### Specific replacements

Providing both the `loc` and `ori` arguments allows specifying the exact configuration used in each replacement.  A zero value for any element of `ori` means to use the optimal replacement at the corresponding location.
"""

# ╔═╡ 183f9ff1-810a-4062-a3c6-cdb82dbd8a7a
begin
    local child =
        substructure_replace(search, replacement; loc=[1, 2, 3, 13], ori=[0, 1, 2, 3])
    view_structure(child)
end

# ╔═╡ 0347c45c-899b-4c4d-9b31-f09a205636f0
md"""
#### Random orientations

By using the `random` keyword argument, the search for optimal alignment can be skipped, and an arbitrary alignment option will be used.
"""

# ╔═╡ b7064b05-6e57-4e81-8cbc-b2075166a1af
begin
    local child = substructure_replace(search, replacement; nb_loc=24, random=true)
    view_structure(child)
end

# ╔═╡ Cell order:
# ╟─359a6c00-c9a6-441d-b258-55bfb5deb4b5
# ╟─8d523993-6e85-443a-9949-12030552b457
# ╠═a948f8b3-4ec5-40b9-b2c1-fcf5b8ad67fa
# ╠═77c63e60-2708-4fef-a6ea-cb394f114b88
# ╟─74b45651-21d5-4332-a4a5-866ea1bb02b8
# ╟─3ff0cd63-8bdf-4ba5-92b6-e2a52f7573c2
# ╠═ccc7fc86-5f48-4cc9-bd5c-349fc1d55b0c
# ╟─b66dd0d5-5bac-4b23-a33f-17305b0d4728
# ╠═a76e79ee-c82c-4771-818d-5380a5fb4c18
# ╟─be86f07e-0669-46c2-8c79-80e65dfcc6f2
# ╠═b53e7c38-d8f5-4f28-a9dc-f7902a86fdb2
# ╟─b402e79a-784b-4d8b-82f1-df4fe1cedad1
# ╟─9b2d534a-4f78-4950-add1-9ba645669bb9
# ╠═4be03110-61ab-4cd6-b3fe-7d51ac5ee771
# ╟─e5eaaef4-13a0-48bd-9f2b-5040b2d10ac1
# ╟─559ef3c3-a176-4d65-8958-810c9b0b32c5
# ╟─d1aa8a19-3702-40f6-b73e-b9ebfc1a7a71
# ╠═44e7e665-ae68-4cd4-b45f-138a0fb8910e
# ╟─af606f02-fbd5-4033-a5f3-56a8f740b1a1
# ╟─c8f8dbd3-4191-460d-944f-f5e456ce8b83
# ╟─127a2bef-3903-4801-bc75-00a6dde2bc6e
# ╟─0e3f7334-9e7f-483b-a0ac-70777902bf51
# ╠═b7b46022-aee0-4d51-8a5e-0c8c005f341a
# ╟─68557426-3204-45a1-8801-84c459e5ac66
# ╠═74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
# ╟─bed9d62f-d61d-4fb1-9a12-188a864fff21
# ╠═9dbf7859-762d-4edc-8300-ac3bba151a8a
# ╟─95215fe9-a7e6-4563-94b6-5ed1d5dde03b
# ╠═8c22f94c-b235-40b2-8fc8-6052c31a6b6e
# ╟─6297833a-ea26-4958-ad56-fc76aeabee69
# ╠═183f9ff1-810a-4062-a3c6-cdb82dbd8a7a
# ╟─0347c45c-899b-4c4d-9b31-f09a205636f0
# ╠═b7064b05-6e57-4e81-8cbc-b2075166a1af
