### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 8ae2593a-cb48-44c9-b9dc-ab901c358a06
begin
    import Pkg
    Pkg.add(url="https://github.com/SimonEnsemble/PoreMatMod.jl")
end

# ╔═╡ 16f0e183-f48d-4dcd-8751-d3c61c875e18
# load required packages (Pluto.jl will automatically install them)
using PoreMatMod, PlutoUI

# ╔═╡ 3cca5289-a829-4717-b796-0834229995d1
begin
    using PoreMatMod.ExampleHelpers
    check_example_data()
end

# ╔═╡ 8d523993-6e85-443a-9949-12030552b457
md"""
## example: generate a hypothetical MOF structure by decorating its linkers with functional groups
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
    we have the IRMOF-1 crystal structure, and wish to append acetamido functional groups to six of its (randomly chosen) BDC (1,4-benzodicarboxylate) linkers to give a mixed-linker IRMOF-1 derivative.

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
    display_query_or_replacement_file("2-!-p-phenylene.xyz")
end

# ╔═╡ d1aa8a19-3702-40f6-b73e-b9ebfc1a7a71
md"""
**Replacement fragment**: next, we construct a replacement fragment as a modified version of the query fragment (with acetamido group in place of one hydrogen atom).
"""

# ╔═╡ 44e7e665-ae68-4cd4-b45f-138a0fb8910e
replacement = moiety("2-acetylamido-p-phenylene.xyz");

# ╔═╡ af606f02-fbd5-4033-a5f3-56a8f740b1a1
view_query_or_replacement("2-acetylamido-p-phenylene.xyz")

# ╔═╡ c8f8dbd3-4191-460d-944f-f5e456ce8b83
with_terminal() do
    display_query_or_replacement_file("2-!-p-phenylene.xyz")
end

# ╔═╡ 127a2bef-3903-4801-bc75-00a6dde2bc6e
md"""
**Find and replace**: finally, we (i) search the parent MOF for the query fragment and (ii) effect the replacements. Ta-da; we have a hypothetical derivatized IRMOF-1 structure. 🎈

n.b. 
* the `nb_loc=6` kwarg indicates that we wish to randomly select 6 matches on the IRMOF-1 parent structure to effect the replacement. the `loc` kwarg grants more control over which BDC linkers are functionalized.
* the acetamido functional groups are appended at random positions on each BDC ligand. the `ori` kwarg grants more control over which positions on each linker are functionalized.
* if we omit `nb_loc=6` as a kwarg, a functional group is appended on all BDC linkers of the parent. 
"""

# ╔═╡ 74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
begin
    child = replace(parent, query => replacement, nb_loc=6)
    view_structure(child)
end

# ╔═╡ 986ecdc4-455f-457e-a964-f00ddfeb53a2
write_cif_message()

# ╔═╡ 71209370-c445-4ff2-a873-b6e31c46419b
write_cif(child, "acetamido-functionalized_IRMOF-1.cif")

# ╔═╡ Cell order:
# ╟─8ae2593a-cb48-44c9-b9dc-ab901c358a06
# ╟─8d523993-6e85-443a-9949-12030552b457
# ╠═16f0e183-f48d-4dcd-8751-d3c61c875e18
# ╠═3cca5289-a829-4717-b796-0834229995d1
# ╠═74b45651-21d5-4332-a4a5-866ea1bb02b8
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
# ╠═74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
# ╟─986ecdc4-455f-457e-a964-f00ddfeb53a2
# ╠═71209370-c445-4ff2-a873-b6e31c46419b
