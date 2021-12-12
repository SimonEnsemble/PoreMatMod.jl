### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 64b72493-6cac-4b92-9a81-6a5df12f7875
begin
    import Pkg
    Pkg.develop("PoreMatMod")
end

# ╔═╡ 5401e009-923e-4a7f-9f3a-fd534f06d8b0
# load required packages (Pluto.jl will automatically install them)
using PoreMatMod, PlutoUI

# ╔═╡ 0f99b225-db96-485c-9e2c-1b4179601e53
using PoreMatMod.ExampleHelpers

# ╔═╡ 8d523993-6e85-443a-9949-12030552b457
md"""
## example: performing a substructure search to find the linkers in a MOF
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
    parent = Crystal("IRMOF-1.cif")    # load .cif file
    infer_bonds!(parent, true)      # infer bonds
    view_structure(parent)          # view structure
end

# ╔═╡ b402e79a-784b-4d8b-82f1-df4fe1cedad1
fragment_construction_note()

# ╔═╡ 9b2d534a-4f78-4950-add1-9ba645669bb9
md"""
**Query fragment**: next, we construct a query fragment as a _p_-phenylene fragment to match that on the BCD linker of the IRMOF-1 parent structure.
"""

# ╔═╡ 4be03110-61ab-4cd6-b3fe-7d51ac5ee771
query = moiety("p-phenylene.xyz");

# ╔═╡ e5eaaef4-13a0-48bd-9f2b-5040b2d10ac1
view_query_or_replacement("p-phenylene.xyz")

# ╔═╡ 559ef3c3-a176-4d65-8958-810c9b0b32c5
with_terminal() do
    display_query_or_replacement_file("p-phenylene.xyz")
end

# ╔═╡ 127a2bef-3903-4801-bc75-00a6dde2bc6e
md"""
### Substructure Search

Find substructures in the `parent` that match the `query`:
"""

# ╔═╡ b7b46022-aee0-4d51-8a5e-0c8c005f341a
search = query ∈ parent

# ╔═╡ 536d3474-bc9c-4576-a120-826020f8d772
hits = isomorphic_substructures(search)

# ╔═╡ ec84ba70-b83e-4a46-9037-aed54bb07b07
view_structure(hits)

# ╔═╡ Cell order:
# ╟─64b72493-6cac-4b92-9a81-6a5df12f7875
# ╟─8d523993-6e85-443a-9949-12030552b457
# ╠═5401e009-923e-4a7f-9f3a-fd534f06d8b0
# ╠═0f99b225-db96-485c-9e2c-1b4179601e53
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
# ╟─127a2bef-3903-4801-bc75-00a6dde2bc6e
# ╠═b7b46022-aee0-4d51-8a5e-0c8c005f341a
# ╠═536d3474-bc9c-4576-a120-826020f8d772
# ╠═ec84ba70-b83e-4a46-9037-aed54bb07b07
