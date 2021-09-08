### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 37939a7a-0651-11ec-11c1-6b5ef0a19ec2
using PoreMatMod, PlutoUI, Bio3DView

# ╔═╡ 21bf7b62-1aef-45c6-9da4-db8d4d69604c
include("ExampleHelper.jl"); # helper functions for viewing molecules and crystals

# ╔═╡ 8d523993-6e85-443a-9949-12030552b457
md"""
## Example: generate a missing-linker defect in UiO-66
"""

# ╔═╡ 8042c752-791d-41a1-b51c-b1073c2e5eee
md"""
The files we read in here are located at:

$(rc[:paths][:crystals])

$(rc[:paths][:moieties])
"""

# ╔═╡ 5b71d14a-be80-4ac3-8983-62571d0d4e7d
md"""
**Task**: We have the crystal structure for UiO-66, and we wish to create a new channel via missing-linker defects with formate ion capping.

**Parent crystal structure**: Below, we read in the `.cif` file describing the parent structure (the UiO-66 unit cell, replicated twice on each crystallographic axis).
"""

# ╔═╡ 8819797c-f1a2-4c46-999d-00316bd21e44
begin
	# read in the parent xtal
	parent = Crystal("UiO-66.cif")	# load .cif file
	infer_bonds!(parent, true)      # infer bonds
	view_structure(parent)          # view structure
end

# ╔═╡ 15ae88e9-b886-4cc6-9ba7-41111bfa06b0
md"""
**Query fragment**: First, we define a query fragment to match what we see in the parent structure.  The atoms to be deleted are tagged: the hydrogens are shown in light pink, and the carbons in dark pink.
"""

# ╔═╡ 1a443428-e283-4205-986e-d0c4ac09bbaa
query = moiety("BDC.xyz");

# ╔═╡ 19282148-e649-4d6e-833d-43fa3cde14c6
view_query_or_replacement("BDC.xyz")

# ╔═╡ d56f8ef4-9a3a-4607-99b1-26034bb23757
with_terminal() do
	display_query_or_replacement_file("BDC.xyz")
end

# ╔═╡ b3a0b763-f9ae-480a-8ad0-ff35f12dc68f
md"""
**Replacement fragment**: Next, we define a replacement fragment as a pair of formate ions, spaced so as to neatly replace the carboxyl groups of the BDC linker.
"""

# ╔═╡ 6a3779f7-dbc5-47b7-a261-1b6144304d5f
replacement = moiety("formate_caps.xyz");

# ╔═╡ 5d62cc3f-4834-4755-b398-922336a26ed8
view_query_or_replacement("formate_caps.xyz")

# ╔═╡ a9005ea7-9fe4-4112-bd9b-a5bb124b0d04
with_terminal() do
	display_query_or_replacement_file("formate_caps.xyz")
end

# ╔═╡ 03b88236-08cb-4f1f-b23d-5f78581373b4
md"""
**Find and replace**: Finally, search the parent MOF for the query fragment and effect the replacements. Nice; we have our UiO-66 with an engineered defect. 🎆
"""

# ╔═╡ 74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
begin
	child = replace(parent, query => replacement; loc=[1, 3, 9])
	view_structure(child)
end

# ╔═╡ Cell order:
# ╟─8d523993-6e85-443a-9949-12030552b457
# ╠═37939a7a-0651-11ec-11c1-6b5ef0a19ec2
# ╠═21bf7b62-1aef-45c6-9da4-db8d4d69604c
# ╟─8042c752-791d-41a1-b51c-b1073c2e5eee
# ╟─5b71d14a-be80-4ac3-8983-62571d0d4e7d
# ╠═8819797c-f1a2-4c46-999d-00316bd21e44
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
