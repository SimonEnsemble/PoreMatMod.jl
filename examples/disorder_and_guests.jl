### A Pluto.jl notebook ###
# v0.16.4

using Markdown
using InteractiveUtils

# ╔═╡ 37939a7a-0651-11ec-11c1-6b5ef0a19ec2
# load required packages (Pluto.jl will automatically install them)
using PoreMatMod, PlutoUI, Bio3DView, Logging

# ╔═╡ 21bf7b62-1aef-45c6-9da4-db8d4d69604c
# helper functions for viewing molecules and crystals and gathering input files
include("ExampleHelper.jl");

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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Bio3DView = "99c8bb3a-9d13-5280-9740-b4880ed9c598"
Logging = "56ddb016-857b-54e1-b83d-db4d58db5568"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PoreMatMod = "2de0d7f0-0963-4438-8bc8-7e7ffe3dc69a"

[compat]
Bio3DView = "~0.1.3"
PlutoUI = "~0.7.14"
PoreMatMod = "~0.2.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BinaryProvider]]
deps = ["Libdl", "Logging", "SHA"]
git-tree-sha1 = "ecdec412a9abc8db54c0efc5548c64dfce072058"
uuid = "b99e7846-7c00-51b0-8f62-c81ae34c0232"
version = "0.5.10"

[[Bio3DView]]
deps = ["Requires"]
git-tree-sha1 = "7f472efd9b6af772307dd017f9deeff2a243754f"
uuid = "99c8bb3a-9d13-5280-9740-b4880ed9c598"
version = "0.1.3"

[[CSV]]
deps = ["Dates", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode"]
git-tree-sha1 = "b83aa3f513be680454437a0eee21001607e5d983"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.8.5"

[[CategoricalArrays]]
deps = ["DataAPI", "Future", "JSON", "Missings", "Printf", "Statistics", "StructTypes", "Unicode"]
git-tree-sha1 = "18d7f3e82c1a80dd38c16453b8fd3f0a7db92f23"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.9.7"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "3533f5a691e60601fe60c90d8bc47a27aa2907ec"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "31d0151f5716b655421d9d75b7fa74cc4e744df2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.39.0"

[[Conda]]
deps = ["JSON", "VersionParsing"]
git-tree-sha1 = "299304989a5e6473d985212c28928899c74e9421"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.5.2"

[[Configurations]]
deps = ["Crayons", "ExproniconLite", "OrderedCollections", "TOML"]
git-tree-sha1 = "b8486a417456d2fbbe2af13e24cef459c9f42429"
uuid = "5218b696-f38b-4ac9-8b61-a12ec717816d"
version = "0.15.4"

[[Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataFrames]]
deps = ["CategoricalArrays", "Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d50972453ef464ddcebdf489d11885468b7b83a3"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "0.22.7"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[ExproniconLite]]
git-tree-sha1 = "c04d5c3442126d75ee4500aa6b0e402cae3bf6ac"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.6.12"

[[FIGlet]]
deps = ["BinaryProvider", "Pkg"]
git-tree-sha1 = "bfc6b52f75b4720581e3e49ae786da6764e65b6a"
uuid = "3064a664-84fe-4d92-92c7-ed492f3d8fae"
version = "0.2.1"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "3c041d2ac0a52a12a27af2782b34900d9c3ee68c"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.11.1"

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FromFile]]
git-tree-sha1 = "f4fce099f9e14cd12d91e29fd105ef891fb55297"
uuid = "ff7dd447-1dcb-4ce3-b8ac-22a812192de7"
version = "0.1.1"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[FuzzyCompletions]]
deps = ["REPL"]
git-tree-sha1 = "2cc2791b324e8ed387a91d7226d17be754e9de61"
uuid = "fb4132e2-a121-4a70-b8a1-d5b831dcdcc2"
version = "0.4.3"

[[GitHubActions]]
deps = ["JSON", "Logging"]
git-tree-sha1 = "56e01ec63d13e1cf015d9ff586156eae3cc7cd6f"
uuid = "6b79fd1a-b13a-48ab-b6b0-aaee1fee41df"
version = "0.1.4"

[[Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "92243c07e786ea3458532e199eb3feee0e7e08eb"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.4.1"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "14eece7a3308b4d8be910e265c724a6ba51a9798"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.16"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "5efcf53d798efede8fee5b2c8b09284be359bf24"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.2"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "f0c6489b12d28fb4c2103073ec7452f3423bd308"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.1"

[[InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLD2]]
deps = ["DataStructures", "FileIO", "MacroTools", "Mmap", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "46b7834ec8165c541b0b5d1c8ba63ec940723ffb"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.15"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "6193c3815f13ba1b78a51ce391db8be016ae9214"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.4"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "5a5bc6bf062f0f95e62d0fe0a2d99699fed82dd9"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.8"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[MetaGraphs]]
deps = ["Graphs", "JLD2", "Random"]
git-tree-sha1 = "5658bebcf7e58e7ff09aa004ff6e806478c5e93a"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.7.0"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f8c673ccc215eb50fcadb285f522420e29e69e1c"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "0.4.5"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "a8cbf066b54d793b9a48c5daa5d586cf2b5bd43d"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.1.0"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "bfd7d8c7fd87f04543810d9cbd3995972236ba1b"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.2"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[Pluto]]
deps = ["Base64", "Configurations", "Dates", "Distributed", "FileWatching", "FuzzyCompletions", "HTTP", "InteractiveUtils", "Logging", "Markdown", "MsgPack", "Pkg", "REPL", "Sockets", "TableIOInterface", "Tables", "UUIDs"]
git-tree-sha1 = "7b4eb6b87276b7ee6604cd9ababf68ba659cfacb"
uuid = "c3e4b0f8-55cb-11ea-2926-15256bba5781"
version = "0.16.4"

[[PlutoSliderServer]]
deps = ["Base64", "Configurations", "Distributed", "FromFile", "GitHubActions", "HTTP", "Logging", "Pkg", "Pluto", "SHA", "Sockets", "TOML", "UUIDs"]
git-tree-sha1 = "cff19495202bd7bfa51fddffca066266913ac05e"
uuid = "2fc8631c-6f24-4c5b-bca7-cbb509c42db4"
version = "0.2.5"

[[PlutoUI]]
deps = ["Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "4c8a7d080daca18545c56f1cac28710c362478f3"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.16"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a193d6ad9c45ada72c14b731a318bedd3c2f00cf"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.3.0"

[[PoreMatMod]]
deps = ["Bio3DView", "CSV", "DataFrames", "FIGlet", "Graphs", "LinearAlgebra", "MetaGraphs", "PlutoSliderServer", "PlutoUI", "Reexport", "StatsBase", "Xtals"]
git-tree-sha1 = "c0e6e7680603c4ad4095a9744893f82301407d2b"
uuid = "2de0d7f0-0963-4438-8bc8-7e7ffe3dc69a"
version = "0.2.4"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "574a6b3ea95f04e8757c0280bb9c29f1a5e35138"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "0.11.1"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "4ba3651d33ef76e24fef6a598b63ffd1c5e1cd17"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.92.5"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "f45b34656397a1f6e729901dc9ef679610bd12b5"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.8"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures", "Random", "Test"]
git-tree-sha1 = "03f5898c9959f8115e30bc7226ada7d0df554ddd"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "0.3.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "eb35dcc66558b2dda84079b9a1be17557d32091a"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.12"

[[StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "d24a825a95a6d98c385001212dc9020d609f2d4f"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.8.1"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableIOInterface]]
git-tree-sha1 = "9a0d3ab8afd14f33a35af7391491ff3104401a35"
uuid = "d1efa939-5518-4425-949f-ab857e148477"
version = "0.1.6"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[VersionParsing]]
git-tree-sha1 = "e575cf85535c7c3292b4d89d89cc29e8c3098e47"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.1"

[[Xtals]]
deps = ["Bio3DView", "CSV", "DataFrames", "Graphs", "JLD2", "LinearAlgebra", "Logging", "MetaGraphs", "Printf", "PyCall", "UUIDs"]
git-tree-sha1 = "51a24da2713c5d8c73e600256d3c4b29df35f067"
uuid = "ede5f01d-793e-4c47-9885-c447d1f18d6d"
version = "0.3.8"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─8d523993-6e85-443a-9949-12030552b457
# ╠═37939a7a-0651-11ec-11c1-6b5ef0a19ec2
# ╠═21bf7b62-1aef-45c6-9da4-db8d4d69604c
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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
