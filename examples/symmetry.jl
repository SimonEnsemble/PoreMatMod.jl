### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ fd82b765-7a62-4f85-bf60-28c83b1731aa
begin
    import Pkg
    Pkg.develop("PoreMatMod")
end

# ╔═╡ 3bb3966e-a7dd-46c1-b649-33169ce424d2
using PoreMatMod, PlutoUI

# ╔═╡ d3eef6f4-ff15-47f3-8686-2c0cb0fb882d
begin
    using PoreMatMod.ExampleHelpers
    check_example_data()
end

# ╔═╡ 8e4fd703-f53e-4056-9466-66f07bacad8d
md"## Example: substructure replacement in a high-symmetry parent"

# ╔═╡ 49ce8a40-3c99-48a8-9403-f5538d1d9a67
input_file_message()

# ╔═╡ 98ee7a1e-1c27-4819-b752-01602cde2799
xtal_folder_message()

# ╔═╡ f5a08f77-9830-4e8e-b746-6a1d61f39009
rc[:paths][:crystals]

# ╔═╡ 8c4c9c65-61a3-4ac8-9413-48aab281237e
moiety_folder_message()

# ╔═╡ 4c5e0230-4db5-49de-b9d7-152ebc834d48
rc[:paths][:moieties]

# ╔═╡ e7be5fba-b120-47f5-a952-8577ee847d99
md"""
!!! example \"the task\"
    We have a crystallographic unit cell of NiPyC. We wish to append methyl groups onto the PyC linkers.

**Parent crystal structure**: first, we read in the `.cif` file describing the parent structure, which is in a high-symmetry representation, meaning this structure must be reproduced according to the symmetry rules of the unit cell's space group. Normally, structures are automatically transformed to P1 symmetry, so we must specify that we wish to preserve the original unit cell.
"""

# ╔═╡ c0e2ffbd-fb85-4187-b8b5-73edea90969a
begin
    # read in the parent xtal, keeping it in its original space group
    parent = Crystal("NiPyC_experiment.cif", convert_to_p1=false) # load cif file
    infer_bonds!(parent, true)                                    # infer bonds
    view_structure(parent)                                        # view structure
end

# ╔═╡ 0bbe1bf4-89e8-4372-8f6a-93fc89334b00
fragment_construction_note()

# ╔═╡ f02fdbf1-2c43-40df-a678-e3c535751f3e
md"""
**Query fragment**: next, we construct a query fragment to match the linker in NiPyC.
"""

# ╔═╡ 16faf8c8-d5c2-4755-945b-546cdac27350
query = moiety("PyC.xyz");

# ╔═╡ 75610c4a-9218-406c-b28f-8e299e197135
view_query_or_replacement("PyC.xyz")

# ╔═╡ 5e90af51-f777-47f9-980b-19bca38de5d4
with_terminal() do
    display_query_or_replacement_file("PyC.xyz")
end

# ╔═╡ 70468997-24fa-4c1b-9f02-379823f97db8
md"""
**Replacement fragment**: then, we construct a replacement fragment as a corrected version of the query fragment (with hydrogen atoms appropriately appended).
"""

# ╔═╡ aac504f5-080e-47b8-9a33-d43c16dc87e7
replacement = moiety("PyC-CH3.xyz");

# ╔═╡ 976795e6-e696-4a37-9cd1-83b4d2695a96
view_query_or_replacement("PyC-CH3.xyz")

# ╔═╡ 7f82be68-558d-4dce-8492-bcc585ddf448
with_terminal() do
    display_query_or_replacement_file("PyC-CH3.xyz")
end

# ╔═╡ 2e3299b5-3b63-48d9-90d7-b8dde5715907
md"""
**Find and replace**: finally, we search the parent MOF for the query fragment and effect the replacement.
"""

# ╔═╡ 0fd8b373-9fec-4b18-b1a6-1b11741e5215
search = query in parent

# ╔═╡ 0b287f09-a121-430a-a749-dc93492d1680
child = substructure_replace(search, replacement, nb_loc=1, wrap=false)

# ╔═╡ d8ea8d5e-f17a-412b-8461-15ba6d9621ec
write_cif(child)

# ╔═╡ eae2225c-40f0-4d68-a9a2-43a39a82f029
view_structure(child)

# ╔═╡ 000144c6-be54-4774-b449-698e9da0741a
write_cif(child, "data/crystals/symmetry_child.cif")

# ╔═╡ 47efaf2e-6758-42c5-aab2-80c1d7725b4a
md"""
**making the super-cell** once we have the high-symmetry structure modified, we can apply the symmetry rules and get a replicated super-cell.
"""

# ╔═╡ 36a1be30-db0e-4c45-8894-859e36793482
begin
    p1_child = Crystal("symmetry_child.cif", check_overlap=false)
    p1_child = replicate(p1_child, (2, 2, 2))
    infer_bonds!(p1_child, true)
    write_cif(p1_child, "supercell.cif")
    view_structure(p1_child)
end

# ╔═╡ Cell order:
# ╟─fd82b765-7a62-4f85-bf60-28c83b1731aa
# ╟─8e4fd703-f53e-4056-9466-66f07bacad8d
# ╠═3bb3966e-a7dd-46c1-b649-33169ce424d2
# ╠═d3eef6f4-ff15-47f3-8686-2c0cb0fb882d
# ╟─49ce8a40-3c99-48a8-9403-f5538d1d9a67
# ╠═98ee7a1e-1c27-4819-b752-01602cde2799
# ╠═f5a08f77-9830-4e8e-b746-6a1d61f39009
# ╟─8c4c9c65-61a3-4ac8-9413-48aab281237e
# ╠═4c5e0230-4db5-49de-b9d7-152ebc834d48
# ╟─e7be5fba-b120-47f5-a952-8577ee847d99
# ╠═c0e2ffbd-fb85-4187-b8b5-73edea90969a
# ╟─0bbe1bf4-89e8-4372-8f6a-93fc89334b00
# ╟─f02fdbf1-2c43-40df-a678-e3c535751f3e
# ╠═16faf8c8-d5c2-4755-945b-546cdac27350
# ╟─75610c4a-9218-406c-b28f-8e299e197135
# ╟─5e90af51-f777-47f9-980b-19bca38de5d4
# ╟─70468997-24fa-4c1b-9f02-379823f97db8
# ╠═aac504f5-080e-47b8-9a33-d43c16dc87e7
# ╟─976795e6-e696-4a37-9cd1-83b4d2695a96
# ╟─7f82be68-558d-4dce-8492-bcc585ddf448
# ╟─2e3299b5-3b63-48d9-90d7-b8dde5715907
# ╠═0fd8b373-9fec-4b18-b1a6-1b11741e5215
# ╠═0b287f09-a121-430a-a749-dc93492d1680
# ╠═d8ea8d5e-f17a-412b-8461-15ba6d9621ec
# ╠═eae2225c-40f0-4d68-a9a2-43a39a82f029
# ╠═000144c6-be54-4774-b449-698e9da0741a
# ╟─47efaf2e-6758-42c5-aab2-80c1d7725b4a
# ╠═36a1be30-db0e-4c45-8894-859e36793482
