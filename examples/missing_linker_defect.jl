### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 53071de7-da2f-4628-9313-2124693ce525

# ╔═╡ d4e77120-8ee6-4514-b066-6127aaa1d6c9
# load required packages (Pluto.jl will automatically install them)
using PoreMatMod, PlutoUI

# ╔═╡ 468c12c0-886f-409b-b36f-b6ff90063e40
begin
    using PoreMatMod.ExampleHelpers
    check_example_data()
end

# ╔═╡ 8d523993-6e85-443a-9949-12030552b457
md"""
## Example: introduce missing-linker defects into a MOF structure
"""

# ╔═╡ e2877131-2e59-4e00-b549-70529d8e71e4
input_file_message()

# ╔═╡ 3eb0313b-0a64-482e-930a-14bd0353a00c
xtal_folder_message()

# ╔═╡ d47c97d5-614e-4c61-b65d-f3fb44014cd1
rc[:paths][:crystals]

# ╔═╡ 528b29d8-847a-4723-8fde-fa319a7b8f3f
moiety_folder_message()

# ╔═╡ 1a62d3bf-95f5-4b1b-af12-6f9db700e5f4
rc[:paths][:moieties]

# ╔═╡ 5b71d14a-be80-4ac3-8983-62571d0d4e7d
md"""
!!! example \"the task\"
    we have the crystal structure for the MOF UiO-66, and we wish to introduce missing-linker defects by removing BDC (1,4-benzodicarboxylate) linkers and adding formate ion capping groups in their place.

**Parent crystal structure**: first, we read in the `.cif` file describing the parent structure (the UiO-66 unit cell, replicated twice on each crystallographic axis).
"""

# ╔═╡ 8819797c-f1a2-4c46-999d-00316bd21e44
begin
    # read in the parent xtal
    parent = Crystal("UiO-66.cif")    # load .cif file
    infer_bonds!(parent, true)      # infer bonds
    view_structure(parent)          # view structure
end

# ╔═╡ 103d3258-5f1c-4b6f-81e3-f45c2bbac844
fragment_construction_note()

# ╔═╡ 15ae88e9-b886-4cc6-9ba7-41111bfa06b0
md"""
**Query fragment**: next, we construct a query fragment to be the BDC linker in the parent structure. the masked atoms to be deleted are tagged with `!` in the `.xyz` input file: the masked hydrogen and carbon atoms are shown in light pink and dark pink, respectively.
"""

# ╔═╡ 1a443428-e283-4205-986e-d0c4ac09bbaa
query = moiety("BDC.xyz");

# ╔═╡ 19282148-e649-4d6e-833d-43fa3cde14c6
view_query_or_replacement("BDC.xyz")

# ╔═╡ d56f8ef4-9a3a-4607-99b1-26034bb23757
with_terminal() do
    return display_query_or_replacement_file("BDC.xyz")
end

# ╔═╡ b3a0b763-f9ae-480a-8ad0-ff35f12dc68f
md"""
**Replacement fragment**: next, we construct a replacement fragment as the BDC linker with the _p_-phenylene ring removed, giving a pair of formate ions, spaced so as to neatly replace (but in effect keep) the carboxyl groups of the BDC linker (and in effect remove the _p_-phenylene ring).
"""

# ╔═╡ 6a3779f7-dbc5-47b7-a261-1b6144304d5f
replacement = moiety("formate_caps.xyz");

# ╔═╡ 5d62cc3f-4834-4755-b398-922336a26ed8
view_query_or_replacement("formate_caps.xyz")

# ╔═╡ a9005ea7-9fe4-4112-bd9b-a5bb124b0d04
with_terminal() do
    return display_query_or_replacement_file("formate_caps.xyz")
end

# ╔═╡ 03b88236-08cb-4f1f-b23d-5f78581373b4
md"""
**Find and replace**: finally, we search the parent MOF for the query fragment and effect the replacements. nice; we have a UiO-66 crystal structure model harboring engineered missing-linker defects. 🎆

n.b.
* the `loc` kwarg indicates which matches to effect the replacement; we carefully chose these locations to introduce a new channel into the MOF.
"""

# ╔═╡ 74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
begin
    child = replace(parent, query => replacement; loc=[3, 9])
    view_structure(child)
end

# ╔═╡ 7b34f300-4815-43f3-8f10-18cb6c76c7f4
write_cif_message()

# ╔═╡ 14f1459b-f505-4950-80f4-7e641abb6b7a
write_cif(child, "defected_UiO-66.cif")

# ╔═╡ a31ddf3a-561c-4240-bb8e-481ed3a5083d
md"""
!!! note \"Regarding fragment selection\"
The astute observer will note that a seemingly simpler pair of query/replacement fragments might be chosen, so as to replace each *half* of a BDC linker with a formate ion.  This is perfectly fine, but leaves the user with the task of determining twice as many locations upon which to operate, such that the correct *six halves* of BDC linkers (making three complete linkers) are replaced.
"""

# ╔═╡ Cell order:
# ╟─53071de7-da2f-4628-9313-2124693ce525
# ╟─8d523993-6e85-443a-9949-12030552b457
# ╠═d4e77120-8ee6-4514-b066-6127aaa1d6c9
# ╠═468c12c0-886f-409b-b36f-b6ff90063e40
# ╟─e2877131-2e59-4e00-b549-70529d8e71e4
# ╠═3eb0313b-0a64-482e-930a-14bd0353a00c
# ╠═d47c97d5-614e-4c61-b65d-f3fb44014cd1
# ╠═528b29d8-847a-4723-8fde-fa319a7b8f3f
# ╠═1a62d3bf-95f5-4b1b-af12-6f9db700e5f4
# ╟─5b71d14a-be80-4ac3-8983-62571d0d4e7d
# ╠═8819797c-f1a2-4c46-999d-00316bd21e44
# ╟─103d3258-5f1c-4b6f-81e3-f45c2bbac844
# ╟─15ae88e9-b886-4cc6-9ba7-41111bfa06b0
# ╠═1a443428-e283-4205-986e-d0c4ac09bbaa
# ╟─19282148-e649-4d6e-833d-43fa3cde14c6
# ╟─d56f8ef4-9a3a-4607-99b1-26034bb23757
# ╟─b3a0b763-f9ae-480a-8ad0-ff35f12dc68f
# ╠═6a3779f7-dbc5-47b7-a261-1b6144304d5f
# ╟─5d62cc3f-4834-4755-b398-922336a26ed8
# ╟─a9005ea7-9fe4-4112-bd9b-a5bb124b0d04
# ╟─03b88236-08cb-4f1f-b23d-5f78581373b4
# ╠═74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
# ╟─7b34f300-4815-43f3-8f10-18cb6c76c7f4
# ╠═14f1459b-f505-4950-80f4-7e641abb6b7a
# ╟─a31ddf3a-561c-4240-bb8e-481ed3a5083d
