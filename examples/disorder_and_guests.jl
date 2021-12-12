### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 8738d1c4-6907-4a7c-96bf-0467b6eb696d
begin
	import Pkg
	Pkg.develop("PoreMatMod")
end

# ╔═╡ ee100c2d-30aa-4b29-b85e-49417e1ed91c
# load required packages (Pluto.jl will automatically install them)
using PoreMatMod, PlutoUI

# ╔═╡ 172f819b-fca8-433a-9a72-e533078e814c
using PoreMatMod.ExampleHelpers

# ╔═╡ 8d523993-6e85-443a-9949-12030552b457
md"""
## Example: correct disorder and remove guest molecules from experimental data
"""

# ╔═╡ 33992496-1700-4033-8deb-694f82bdd1bc
input_file_message()

# ╔═╡ dc03627f-3da3-4175-b52e-664522df2302
xtal_folder_message()

# ╔═╡ 4be347c1-3f1a-40c1-ae38-f82ed692381e
rc[:paths][:crystals]

# ╔═╡ de84b4a3-a365-4860-9729-2377367c64db
moiety_folder_message()

# ╔═╡ b6dff24a-d0a0-4824-8bb2-2dcb60196b0a
rc[:paths][:moieties]

# ╔═╡ 5b71d14a-be80-4ac3-8983-62571d0d4e7d
md"""
!!! example \"the task\"
	we have the experimental structure of SIFSIX-Cu-2-i, which features disordered pyridyl rings in the linkers---presumably as an artifact of X-ray structure determination---and acetylene guest molecules in the pores.  We wish to (i) correct the disorder by selecting a single conformation for each ring and (ii) remove the guest molecules from its pores.

**Parent crystal structure**: first, we read in the `.cif` file describing the SIFSIX-Cu-2-i parent structure, which presents disordered ligands and guest molecules in the pores.
"""

# ╔═╡ 6d9c3b97-fd26-470c-8acf-8ce10b33b82d
begin
	# read in the parent xtal
	parent = Crystal("SIFSIX-2-Cu-i.cif", check_overlap=false) 
	
	infer_bonds!(parent, true)	# infer bonds
	view_structure(parent)      # view structure
end

# ╔═╡ c0feae1a-a228-40eb-a234-e58617fa71dd
fragment_construction_note()

# ╔═╡ ba265b61-9341-48fd-9f93-befd3e65d315
md"""
**Query fragment 1 (disordered ring)**: second, we construct a query fragment to match the disordered rings in the parent structure by cutting one out of the parent structure. we mask all atoms of the query fragment except those belonging to a single conformation of the ring. the masked hydrogen and carbon atoms below are shown in light pink and dark pink, respectively.
"""

# ╔═╡ f2907335-798c-4f9f-bbab-fa8763fb43bb
query_disordered_ring = moiety("disordered_ligand!.xyz");

# ╔═╡ efd5d374-929c-4851-9219-d8e2b86ebb85
view_query_or_replacement("disordered_ligand!.xyz")

# ╔═╡ 3077d0db-01e7-4aab-b183-c0f69a1f2da3
with_terminal() do
	display_query_or_replacement_file("disordered_ligand!.xyz")
end

# ╔═╡ 10911157-64eb-485b-86f1-e83dc201b054
md"""
**Query fragment 2 (guest molecule)**: we also construct an acetylene query fragment to match the guest molecules in the parent structure.
"""

# ╔═╡ db0efdc0-9d29-46ca-8992-39cd7a9ad36c
query_guest = moiety("acetylene.xyz");

# ╔═╡ b5a090eb-5796-44bf-a19c-3705e99d26dd
view_query_or_replacement("acetylene.xyz")

# ╔═╡ 83d00463-ae3b-47d8-b65d-ff8a0aa522e8
with_terminal() do
	display_query_or_replacement_file("acetylene.xyz")
end

# ╔═╡ 25ed4da3-525c-4045-82c4-b3cbf8e707e3
md"""
**Replacement fragment**: next, we construct a replacement fragment---a (non-disordered) pyridyl ring---as a corrected version of the first query fragment that represented the disordered ligand.
"""

# ╔═╡ b57e57d8-897b-466a-b2b5-517afd969123
replacement_ring = moiety("4-pyridyl.xyz");

# ╔═╡ fee26c78-75c2-4844-acdb-d2c4cd71d654
view_query_or_replacement("4-pyridyl.xyz")

# ╔═╡ e85d27a8-ec86-4e38-92a4-f8a2ad9a0ce3
with_terminal() do
	display_query_or_replacement_file("4-pyridyl.xyz")
end

# ╔═╡ da89adc2-0688-44c6-9fd2-791bd13c8d74
md"""
**Find and replace**: finally, 
1. search the parent MOF for the first query fragment (disordered rings) and effect the replacements with the corrected (single-conformation ring), then
2. search the parent for the second query fragment (guest molecules) as disconnected components, and delete them.

Hooray; we have a simulation-ready SIFSIX-Cu-2-i structure. 💖

n.b.
* if we had not passed the `disconnected_component=true` kwarg, the algorithm would have found and removed the sides of the pyridyl rings which contain a H-C-C-H fragment!
"""

# ╔═╡ 74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
begin
	# repair ring disorder
	child = replace(parent, query_disordered_ring => replacement_ring)
	
	# search for disconnected acetylene components
	search = substructure_search(query_guest, child, disconnected_component=true)
	# delete guest molecules
	child = substructure_replace(search, nothing)
	
	view_structure(child)
end

# ╔═╡ 11095805-7c76-42f4-836d-919c2cb27d1c
write_cif_message()

# ╔═╡ c5e14b94-427d-4684-9dd7-bc44d965c570
write_cif(child, "simulation_ready_SIFSIX-Cu-2-i.cif")

# ╔═╡ Cell order:
# ╟─8738d1c4-6907-4a7c-96bf-0467b6eb696d
# ╟─8d523993-6e85-443a-9949-12030552b457
# ╠═ee100c2d-30aa-4b29-b85e-49417e1ed91c
# ╠═172f819b-fca8-433a-9a72-e533078e814c
# ╟─33992496-1700-4033-8deb-694f82bdd1bc
# ╟─dc03627f-3da3-4175-b52e-664522df2302
# ╠═4be347c1-3f1a-40c1-ae38-f82ed692381e
# ╟─de84b4a3-a365-4860-9729-2377367c64db
# ╠═b6dff24a-d0a0-4824-8bb2-2dcb60196b0a
# ╟─5b71d14a-be80-4ac3-8983-62571d0d4e7d
# ╠═6d9c3b97-fd26-470c-8acf-8ce10b33b82d
# ╟─c0feae1a-a228-40eb-a234-e58617fa71dd
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
# ╟─11095805-7c76-42f4-836d-919c2cb27d1c
# ╠═c5e14b94-427d-4684-9dd7-bc44d965c570
