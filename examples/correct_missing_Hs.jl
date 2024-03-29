### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 0069ec00-f42b-40bb-a82e-c91ec78583e4

# ╔═╡ 37939a7a-0651-11ec-11c1-6b5ef0a19ec2
# load required packages (Pluto.jl will automatically install them)
using PoreMatMod, PlutoUI

# ╔═╡ 2889e30d-442d-4620-b5f6-0216e21c623b
begin
    using PoreMatMod.ExampleHelpers
    check_example_data()
end

# ╔═╡ 8d523993-6e85-443a-9949-12030552b457
md"""
## Example: append missing hydrogen atoms to linkers of a MOF
"""

# ╔═╡ ab9b556e-5082-48fa-a7ca-b37ef895d703
input_file_message()

# ╔═╡ 9a9bd9fd-5fea-4c01-a87a-5a5a0d332c1a
xtal_folder_message()

# ╔═╡ fcebd06c-9607-42a5-ba1b-e7683b0924ab
rc[:paths][:crystals]

# ╔═╡ af6f9a09-17ad-45dd-b1fd-50f9cc5a692e
moiety_folder_message()

# ╔═╡ 62439d74-601c-40fc-a488-c0838b9ada69
rc[:paths][:moieties]

# ╔═╡ 026100b5-0708-48bb-840d-931605524874
md"""
!!! example \"the task\"
    We have an IRMOF-1 crystal structure with hydrogen atoms missing on the linkers, presumably owing to artifacts of X-ray structure determination. We wish to append hydrogen atoms onto the missing positions on the linkers.

**Parent crystal structure**: first, we read in the `.cif` file describing the parent structure, which is not simulation-ready owing to the missing hydrogen atoms.
"""

# ╔═╡ 0433da26-4f59-424f-9603-875d904c0fd5
begin
    # read in the parent xtal
    parent = Crystal("IRMOF-1_noH.cif") # load .cif file
    infer_bonds!(parent, true)          # infer bonds
    view_structure(parent)              # view structure
end

# ╔═╡ 4ad53cf5-d063-4bb6-ae20-1b6cc29902c8
fragment_construction_note()

# ╔═╡ 64b95411-07da-44af-a06e-9e6676328ffd
md"""
**Query fragment**: next, we construct a query fragment to match the corrupted fragment in the IRMOF-1 parent structure.
"""

# ╔═╡ 83532414-6471-4002-b23c-1600243318d1
query = moiety("1,4-C-phenylene_noH.xyz");

# ╔═╡ fc9e8e21-02c0-43ca-980f-55496526d7f3
view_query_or_replacement("1,4-C-phenylene_noH.xyz")

# ╔═╡ f68335d7-b4e0-40b3-b10d-bf406ab42c1c
with_terminal() do
    return display_query_or_replacement_file("1,4-C-phenylene_noH.xyz")
end

# ╔═╡ c53f896d-fa27-4290-aa6d-aa8c0c467f3b
md"""
**Replacement fragment**: then, we construct a replacement fragment as a corrected version of the query fragment (with hydrogen atoms appropriately appended).
"""

# ╔═╡ 0d8ac0fd-c7b1-4781-aaa2-33cc8c1c08ae
replacement = moiety("1,4-C-phenylene.xyz");

# ╔═╡ 836f6d08-c507-44c6-b927-c9ea240f40f8
view_query_or_replacement("1,4-C-phenylene.xyz")

# ╔═╡ b172c36b-80bf-4620-b2e4-5c39d719962e
with_terminal() do
    return display_query_or_replacement_file("1,4-C-phenylene.xyz")
end

# ╔═╡ 5b71d14a-be80-4ac3-8983-62571d0d4e7d
md"""
**Find and replace**: finally, we search the parent MOF for the query fragment and effect the replacements. Voila; we have a simulation-ready IRMOF-1 structure. 🚀
"""

# ╔═╡ 74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
begin
    child = replace(parent, query => replacement)
    view_structure(child)
end

# ╔═╡ 5869bf7a-958e-4e51-997a-18497e7deaba
write_cif_message()

# ╔═╡ dbf3a196-dd4d-4ac8-8396-42ac7c7e0ba1
write_cif(child, "simulation_ready_IRMOF-1.cif")

# ╔═╡ Cell order:
# ╠═0069ec00-f42b-40bb-a82e-c91ec78583e4
# ╟─8d523993-6e85-443a-9949-12030552b457
# ╠═37939a7a-0651-11ec-11c1-6b5ef0a19ec2
# ╠═2889e30d-442d-4620-b5f6-0216e21c623b
# ╟─ab9b556e-5082-48fa-a7ca-b37ef895d703
# ╟─9a9bd9fd-5fea-4c01-a87a-5a5a0d332c1a
# ╠═fcebd06c-9607-42a5-ba1b-e7683b0924ab
# ╟─af6f9a09-17ad-45dd-b1fd-50f9cc5a692e
# ╠═62439d74-601c-40fc-a488-c0838b9ada69
# ╟─026100b5-0708-48bb-840d-931605524874
# ╠═0433da26-4f59-424f-9603-875d904c0fd5
# ╟─4ad53cf5-d063-4bb6-ae20-1b6cc29902c8
# ╟─64b95411-07da-44af-a06e-9e6676328ffd
# ╠═83532414-6471-4002-b23c-1600243318d1
# ╟─fc9e8e21-02c0-43ca-980f-55496526d7f3
# ╟─f68335d7-b4e0-40b3-b10d-bf406ab42c1c
# ╟─c53f896d-fa27-4290-aa6d-aa8c0c467f3b
# ╠═0d8ac0fd-c7b1-4781-aaa2-33cc8c1c08ae
# ╟─836f6d08-c507-44c6-b927-c9ea240f40f8
# ╟─b172c36b-80bf-4620-b2e4-5c39d719962e
# ╟─5b71d14a-be80-4ac3-8983-62571d0d4e7d
# ╠═74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
# ╟─5869bf7a-958e-4e51-997a-18497e7deaba
# ╠═dbf3a196-dd4d-4ac8-8396-42ac7c7e0ba1
