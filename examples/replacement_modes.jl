### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# ╔═╡ 37939a7a-0651-11ec-11c1-6b5ef0a19ec2
# load required packages (Pluto.jl will automatically install them)
using PoreMatMod, PlutoUI, Bio3DView, Logging

# ╔═╡ 8ca1eb06-df90-4837-87b1-2e76e1670504
# helper functions for viewing molecules and crystals and gathering input files
include("ExampleHelper.jl");

# ╔═╡ 8d523993-6e85-443a-9949-12030552b457
md"""
## example: generating hypothetical MOF structures with different replacement styles
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
	parent = Crystal("IRMOF-1.cif")	# load .cif file
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
replacement = moiety("2-nitro-p-phenylene.xyz");

# ╔═╡ af606f02-fbd5-4033-a5f3-56a8f740b1a1
view_query_or_replacement("2-nitro-p-phenylene.xyz")

# ╔═╡ c8f8dbd3-4191-460d-944f-f5e456ce8b83
with_terminal() do
	display_query_or_replacement_file("2-nitro-p-phenylene.xyz")
end

# ╔═╡ 127a2bef-3903-4801-bc75-00a6dde2bc6e
md"""
### Replacement Modes

There are several options when performing a replacement that control how the replacement fragments will be positioned in the parent structure.

With all three file inputs loaded (IRMOF-1 as `parent`, 2-!-*p*-phenylene as `query`, and 2-nitro-*p*-phenylene as `replacement`) and a `search` performed, replacements can be made.

`PoreMatMod.jl` has several replacement modes, one of which must be specified.
"""

# ╔═╡ 0e3f7334-9e7f-483b-a0ac-70777902bf51
md"""
!!! note \"Note\"
	Multiple replacements can be done with a single search.
"""

# ╔═╡ b7b46022-aee0-4d51-8a5e-0c8c005f341a
search = query ∈ parent

# ╔═╡ 68557426-3204-45a1-8801-84c459e5ac66
md"""
#### Default: optimal orientation at all locations

Optimal configurations will be chosen for each location in `search.results`, so that each occurrence of the `query` in the `parent` is replaced with minimal perturbation of the conserved atoms from the parent structure.
"""

# ╔═╡ 74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
begin
	local child = substructure_replace(search, replacement)
	view_structure(child)
end

# ╔═╡ bed9d62f-d61d-4fb1-9a12-188a864fff21
md"""
#### Optimal Orientation at *n* random locations

If only some of the linkers should be replaced, the `nb_loc` argument lets us specify how many.
"""

# ╔═╡ 9dbf7859-762d-4edc-8300-ac3bba151a8a
begin
	local child = substructure_replace(search, replacement, nb_loc=8)
	view_structure(child)
end

# ╔═╡ 95215fe9-a7e6-4563-94b6-5ed1d5dde03b
md"""
#### Optimal orientation at specific locations

Specific locations are chosen by providing the `loc` argument.
"""

# ╔═╡ 8c22f94c-b235-40b2-8fc8-6052c31a6b6e
begin
	local child = substructure_replace(search, replacement, loc=[17, 18, 19, 20])
	view_structure(child)
end

# ╔═╡ 6297833a-ea26-4958-ad56-fc76aeabee69
md"""
#### Specific replacements

Providing both the `loc` and `ori` arguments allows specifying the exact configuration used in each replacement.  A zero value for any element of `ori` means to use the optimal replacement at the corresponding location.
"""

# ╔═╡ 183f9ff1-810a-4062-a3c6-cdb82dbd8a7a
begin
	local child = substructure_replace(
		search, replacement, loc=[1, 2, 3, 13], ori=[0, 1, 2, 3])
	view_structure(child)
end

# ╔═╡ 0347c45c-899b-4c4d-9b31-f09a205636f0
md"""
#### Random orientations

By using the `random` keyword argument, the search for optimal alignment can be skipped, and an arbitrary alignment option will be used.
"""

# ╔═╡ b7064b05-6e57-4e81-8cbc-b2075166a1af
begin
	local child = substructure_replace(search, replacement, nb_loc=24, random=true)
	view_structure(child)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Bio3DView = "99c8bb3a-9d13-5280-9740-b4880ed9c598"
Logging = "56ddb016-857b-54e1-b83d-db4d58db5568"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PoreMatMod = "2de0d7f0-0963-4438-8bc8-7e7ffe3dc69a"

[compat]
Bio3DView = "~0.1.3"
PlutoUI = "~0.7.12"
PoreMatMod = "~0.1.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AssetRegistry]]
deps = ["Distributed", "JSON", "Pidfile", "SHA", "Test"]
git-tree-sha1 = "b25e88db7944f98789130d7b503276bc34bc098e"
uuid = "bf4720bc-e11a-5d0c-854e-bdca1663c893"
version = "0.1.0"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BinDeps]]
deps = ["Libdl", "Pkg", "SHA", "URIParser", "Unicode"]
git-tree-sha1 = "1289b57e8cf019aede076edab0587eb9644175bd"
uuid = "9e28174c-4ba2-5203-b857-d8d62c4213ee"
version = "1.0.2"

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

[[Blink]]
deps = ["Base64", "BinDeps", "Distributed", "JSExpr", "JSON", "Lazy", "Logging", "MacroTools", "Mustache", "Mux", "Reexport", "Sockets", "WebIO", "WebSockets"]
git-tree-sha1 = "08d0b679fd7caa49e2bca9214b131289e19808c0"
uuid = "ad839575-38b3-5650-b840-f874b8c74a25"
version = "0.12.5"

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

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

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

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[ExproniconLite]]
git-tree-sha1 = "45b421f664eba4cd4bed8e5e706f69f2ccd8fcb3"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.6.10"

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

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FromFile]]
git-tree-sha1 = "f4fce099f9e14cd12d91e29fd105ef891fb55297"
uuid = "ff7dd447-1dcb-4ce3-b8ac-22a812192de7"
version = "0.1.1"

[[FunctionalCollections]]
deps = ["Test"]
git-tree-sha1 = "04cb9cfaa6ba5311973994fe3496ddec19b6292a"
uuid = "de31a74c-ac4f-5751-b3fd-e18cd04993ca"
version = "0.5.0"

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

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "24675428ca27678f003414a98c9e473e45fe6a21"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.15"

[[Hiccup]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "6187bb2d5fcbb2007c39e7ac53308b0d371124bd"
uuid = "9fb69e20-1954-56bb-a84f-559cc56a8ff7"
version = "0.2.2"

[[HypertextLiteral]]
git-tree-sha1 = "72053798e1be56026b81d4e2682dbe58922e5ec9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.0"

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

[[InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLD2]]
deps = ["DataStructures", "FileIO", "MacroTools", "Mmap", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "192934b3e2a94e897ce177423fd6cf7bdf464bce"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.14"

[[JSExpr]]
deps = ["JSON", "MacroTools", "Observables", "WebIO"]
git-tree-sha1 = "bd6c034156b1e7295450a219c4340e32e50b08b1"
uuid = "97c1335a-c9c5-57fe-bc5d-ec35cebe8660"
version = "0.5.3"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

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

[[LightGraphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "432428df5f360964040ed60418dd5601ecd240b6"
uuid = "093fc24a-ae57-5d10-9952-331d41423f4d"
version = "1.3.5"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

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
deps = ["JLD2", "LightGraphs", "Random"]
git-tree-sha1 = "81c0488104fb0dc977f38b4acaff81e6e79efc4d"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.6.8"

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

[[Mustache]]
deps = ["Printf", "Tables"]
git-tree-sha1 = "36995ef0d532fe08119d70b2365b7b03d4e00f48"
uuid = "ffc61752-8dc7-55ee-8c37-f3e9cdd09e70"
version = "1.0.10"

[[Mux]]
deps = ["AssetRegistry", "Base64", "HTTP", "Hiccup", "Pkg", "Sockets", "WebSockets"]
git-tree-sha1 = "82dfb2cead9895e10ee1b0ca37a01088456c4364"
uuid = "a975b10e-0019-58db-a62f-e48ff68538c9"
version = "0.7.6"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Observables]]
git-tree-sha1 = "fe29afdef3d0c4a8286128d4e45cc50621b1e43d"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.4.0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "bfd7d8c7fd87f04543810d9cbd3995972236ba1b"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.2"

[[Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "1be8660b2064893cd2dae4bd004b589278e4440d"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.2.0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[Pluto]]
deps = ["Base64", "Configurations", "Dates", "Distributed", "FileWatching", "FuzzyCompletions", "HTTP", "InteractiveUtils", "Logging", "Markdown", "MsgPack", "Pkg", "REPL", "Sockets", "TableIOInterface", "Tables", "UUIDs"]
git-tree-sha1 = "d4955d6b5267ed826d03e3f05cc23426492b23c9"
uuid = "c3e4b0f8-55cb-11ea-2926-15256bba5781"
version = "0.16.1"

[[PlutoSliderServer]]
deps = ["Base64", "Configurations", "Distributed", "FromFile", "GitHubActions", "HTTP", "Logging", "Pkg", "Pluto", "SHA", "Sockets", "TOML", "UUIDs"]
git-tree-sha1 = "cff19495202bd7bfa51fddffca066266913ac05e"
uuid = "2fc8631c-6f24-4c5b-bca7-cbb509c42db4"
version = "0.2.5"

[[PlutoUI]]
deps = ["Base64", "Dates", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "f35ae11e070dbf123d5a6f54cbda45818d765ad2"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.12"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a193d6ad9c45ada72c14b731a318bedd3c2f00cf"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.3.0"

[[PoreMatMod]]
deps = ["Bio3DView", "Blink", "CSV", "DataFrames", "FIGlet", "LightGraphs", "LinearAlgebra", "MetaGraphs", "PlutoSliderServer", "PlutoUI", "Reexport", "StatsBase", "Xtals"]
git-tree-sha1 = "4234ac253ba7c16a30ee8dc562960e509ff456c6"
uuid = "2de0d7f0-0963-4438-8bc8-7e7ffe3dc69a"
version = "0.1.2"

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
git-tree-sha1 = "169bb8ea6b1b143c5cf57df6d34d022a7b60c6db"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.92.3"

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
git-tree-sha1 = "54f37736d8934a12a200edea2f9206b03bdf3159"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.7"

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
git-tree-sha1 = "3240808c6d463ac46f1c1cd7638375cd22abbccb"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.12"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8cbbc098554648c84f79a463c9ff0fd277144b6c"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.10"

[[StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "8445bf99a36d703a09c601f9a57e2f83000ef2ae"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.7.3"

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
git-tree-sha1 = "1162ce4a6c4b7e31e0e6b14486a6986951c73be9"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.5.2"

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

[[URIParser]]
deps = ["Unicode"]
git-tree-sha1 = "53a9f49546b8d2dd2e688d216421d050c9a31d0d"
uuid = "30578b45-9adc-5946-b283-645ec420af67"
version = "0.4.1"

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
git-tree-sha1 = "80229be1f670524750d905f8fc8148e5a8c4537f"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.0"

[[WebIO]]
deps = ["AssetRegistry", "Base64", "Distributed", "FunctionalCollections", "JSON", "Logging", "Observables", "Pkg", "Random", "Requires", "Sockets", "UUIDs", "WebSockets", "Widgets"]
git-tree-sha1 = "5fe32e4086d49f7ab9b087296742859f3ae6d62a"
uuid = "0f1e0344-ec1d-5b48-a673-e5cf874b6c29"
version = "0.8.16"

[[WebSockets]]
deps = ["Base64", "Dates", "HTTP", "Logging", "Sockets"]
git-tree-sha1 = "f91a602e25fe6b89afc93cf02a4ae18ee9384ce3"
uuid = "104b5d7c-a370-577a-8038-80a2059c5097"
version = "1.5.9"

[[Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "80661f59d28714632132c73779f8becc19a113f2"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.4"

[[Xtals]]
deps = ["Bio3DView", "CSV", "DataFrames", "JLD2", "LightGraphs", "LinearAlgebra", "Logging", "MetaGraphs", "Printf", "PyCall", "UUIDs"]
git-tree-sha1 = "c2cf3eea05391f316210e3b9b6037e5b4721f8c6"
uuid = "ede5f01d-793e-4c47-9885-c447d1f18d6d"
version = "0.3.4"

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
# ╠═8ca1eb06-df90-4837-87b1-2e76e1670504
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
# ╟─d1aa8a19-3702-40f6-b73e-b9ebfc1a7a71
# ╠═44e7e665-ae68-4cd4-b45f-138a0fb8910e
# ╟─af606f02-fbd5-4033-a5f3-56a8f740b1a1
# ╟─c8f8dbd3-4191-460d-944f-f5e456ce8b83
# ╟─127a2bef-3903-4801-bc75-00a6dde2bc6e
# ╟─0e3f7334-9e7f-483b-a0ac-70777902bf51
# ╠═b7b46022-aee0-4d51-8a5e-0c8c005f341a
# ╟─68557426-3204-45a1-8801-84c459e5ac66
# ╠═74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
# ╟─bed9d62f-d61d-4fb1-9a12-188a864fff21
# ╠═9dbf7859-762d-4edc-8300-ac3bba151a8a
# ╟─95215fe9-a7e6-4563-94b6-5ed1d5dde03b
# ╠═8c22f94c-b235-40b2-8fc8-6052c31a6b6e
# ╟─6297833a-ea26-4958-ad56-fc76aeabee69
# ╠═183f9ff1-810a-4062-a3c6-cdb82dbd8a7a
# ╟─0347c45c-899b-4c4d-9b31-f09a205636f0
# ╠═b7064b05-6e57-4e81-8cbc-b2075166a1af
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
