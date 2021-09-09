### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 37939a7a-0651-11ec-11c1-6b5ef0a19ec2
begin
	import Pkg
	Pkg.activate("/home/cokes/.julia/dev/PoreMatMod")
	using PoreMatMod, PlutoUI, Bio3DView
end

# ╔═╡ 8ca1eb06-df90-4837-87b1-2e76e1670504
include("ExampleHelper.jl"); # helper functions for viewing molecules and crystals

# ╔═╡ 8d523993-6e85-443a-9949-12030552b457
md"""
## Example: generate a hypothetical MOF structure by replacing some of the linkers in IRMOF-1
"""

# ╔═╡ 5b71d14a-be80-4ac3-8983-62571d0d4e7d
md"""
The files we read in here are located at:

$(rc[:paths][:crystals])

$(rc[:paths][:moieties])
"""

# ╔═╡ be86f07e-0669-46c2-8c79-80e65dfcc6f2
md"""
**Task**: We have the IRMOF-1 crystal structure, and wish to replace a certain proportion of the linkers (BDC) with a derivative (acetamido-BDC).

**Parent crystal structure**: Below, we read in the `.cif` file describing the parent structure.
"""

# ╔═╡ b53e7c38-d8f5-4f28-a9dc-f7902a86fdb2
begin
	# read in the parent xtal
	parent = Crystal("IRMOF-1.cif")	# load .cif file
	infer_bonds!(parent, true)      # infer bonds
	view_structure(parent)          # view structure
end

# ╔═╡ 9b2d534a-4f78-4950-add1-9ba645669bb9
md"""
**Query fragment**: First, we define a query fragment to match what we see in the parent structure, with one hydrogen atom "masked" by tagging its species label with `!` in the input. The masked hydrogen is shown in pink.
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
**Replacement fragment**: Next, we define a replacement fragment as a modified version of the query fragment (with acetamido group in place of one hydrogen).
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
**Find and replace**: Finally, search the parent MOF for the query fragment and effect the replacements. Ta-da; we have a hypothetical derivatized IRMOF-1 structure. 🎈
"""

# ╔═╡ 74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
begin
	child = replace(parent, query => replacement, nb_loc=6)
	view_structure(child)
end

# ╔═╡ Cell order:
# ╟─8d523993-6e85-443a-9949-12030552b457
# ╠═37939a7a-0651-11ec-11c1-6b5ef0a19ec2
# ╠═8ca1eb06-df90-4837-87b1-2e76e1670504
# ╟─5b71d14a-be80-4ac3-8983-62571d0d4e7d
# ╟─be86f07e-0669-46c2-8c79-80e65dfcc6f2
# ╠═b53e7c38-d8f5-4f28-a9dc-f7902a86fdb2
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
