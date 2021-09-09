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
## Example: correct disorder and remove guest molecules from experimental data
"""

# ╔═╡ 7d08c908-87dc-4e57-85d6-74ba923a33e4
md"""
The files we read in here are located at:

$(rc[:paths][:crystals])

$(rc[:paths][:moieties])
"""

# ╔═╡ 5b71d14a-be80-4ac3-8983-62571d0d4e7d
md"""
**Task**: We have the experimental structure of SIFSIX-Cu-2-i, which features disordered pyridyl rings in the linkers and acetylene guest molecules in the pores.  We wish to remove the disorder by selecting a single state for each ring, and delete the guest molecules.

**Parent crystal structure**: Below, we read in the `.cif` file describing the experimental parent structure.
"""

# ╔═╡ 6d9c3b97-fd26-470c-8acf-8ce10b33b82d
begin
	# read in the parent xtal
	
	# load .cif file
	parent = Crystal("EMEHUB_C2H2.cif", remove_duplicates=true, check_overlap=false) 
	
	infer_bonds!(parent, true)	# infer bonds
	view_structure(parent)      # view structure
end

# ╔═╡ ba265b61-9341-48fd-9f93-befd3e65d315
md"""
**Query fragment 1**: First, we define a query fragment to match what we see for the disordered rings in the parent structure, with only a single state of the ring un-masked. Masked hydrogen atoms are shown in light pink; masked carbons in dark pink.
"""

# ╔═╡ f2907335-798c-4f9f-bbab-fa8763fb43bb
query1 = moiety("disordered_ligand!.xyz");

# ╔═╡ efd5d374-929c-4851-9219-d8e2b86ebb85
view_query_or_replacement("disordered_ligand!.xyz")

# ╔═╡ 3077d0db-01e7-4aab-b183-c0f69a1f2da3
with_terminal() do
	display_query_or_replacement_file("disordered_ligand!.xyz")
end

# ╔═╡ 10911157-64eb-485b-86f1-e83dc201b054
md"""
**Query fragment 2**: We must also define a query fragment to match what we see for the guest molecules in the parent structure.
"""

# ╔═╡ db0efdc0-9d29-46ca-8992-39cd7a9ad36c
query2 = moiety("acetylene.xyz");

# ╔═╡ b5a090eb-5796-44bf-a19c-3705e99d26dd
view_query_or_replacement("acetylene.xyz")

# ╔═╡ 83d00463-ae3b-47d8-b65d-ff8a0aa522e8
with_terminal() do
	display_query_or_replacement_file("acetylene.xyz")
end

# ╔═╡ 25ed4da3-525c-4045-82c4-b3cbf8e707e3
md"""
**Replacement fragment**: Next, we define a replacement fragment as a corrected version of the first query fragment (a pyridyl ring).
"""

# ╔═╡ b57e57d8-897b-466a-b2b5-517afd969123
replacement = moiety("4-pyridyl.xyz");

# ╔═╡ fee26c78-75c2-4844-acdb-d2c4cd71d654
view_query_or_replacement("4-pyridyl.xyz")

# ╔═╡ e85d27a8-ec86-4e38-92a4-f8a2ad9a0ce3
with_terminal() do
	display_query_or_replacement_file("4-pyridyl.xyz")
end

# ╔═╡ da89adc2-0688-44c6-9fd2-791bd13c8d74
md"""
**Find and replace**: Finally, search the parent MOF for the first query fragment (disordered rings) and effect the replacements; then search the parent for the second query fragment (guest molecules) as disconnected components, and delete them. Hooray; we have a simulation-ready SIFSIX-Cu-2-i structure. 💖
"""

# ╔═╡ 74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
begin
	# repair ring disorder
	child = replace(parent, query1 => replacement)
	
	# search for disconnected acetylene components
	search = substructure_search(query2, child, disconnected_component=true)
	# delete guest molecules
	child = substructure_replace(search, nothing)
	
	view_structure(child)
end

# ╔═╡ Cell order:
# ╟─8d523993-6e85-443a-9949-12030552b457
# ╠═37939a7a-0651-11ec-11c1-6b5ef0a19ec2
# ╠═21bf7b62-1aef-45c6-9da4-db8d4d69604c
# ╟─7d08c908-87dc-4e57-85d6-74ba923a33e4
# ╟─5b71d14a-be80-4ac3-8983-62571d0d4e7d
# ╠═6d9c3b97-fd26-470c-8acf-8ce10b33b82d
# ╟─ba265b61-9341-48fd-9f93-befd3e65d315
# ╠═f2907335-798c-4f9f-bbab-fa8763fb43bb
# ╟─efd5d374-929c-4851-9219-d8e2b86ebb85
# ╟─3077d0db-01e7-4aab-b183-c0f69a1f2da3
# ╟─10911157-64eb-485b-86f1-e83dc201b054
# ╠═db0efdc0-9d29-46ca-8992-39cd7a9ad36c
# ╟─b5a090eb-5796-44bf-a19c-3705e99d26dd
# ╟─83d00463-ae3b-47d8-b65d-ff8a0aa522e8
# ╟─25ed4da3-525c-4045-82c4-b3cbf8e707e3
# ╠═b57e57d8-897b-466a-b2b5-517afd969123
# ╟─fee26c78-75c2-4844-acdb-d2c4cd71d654
# ╟─e85d27a8-ec86-4e38-92a4-f8a2ad9a0ce3
# ╟─da89adc2-0688-44c6-9fd2-791bd13c8d74
# ╠═74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
