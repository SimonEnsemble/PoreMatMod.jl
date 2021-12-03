### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 9c36b2ed-28f1-492c-949d-4bf0e871152b
begin
	import Pkg
	Pkg.develop("PoreMatMod")
end

# ╔═╡ 42005e40-5450-11ec-3384-69ab8820daa2
using PoreMatMod, Bio3DView

# ╔═╡ 2b10d534-0c82-4c5f-8233-82ab7c519ce1
include("ExampleHelper.jl")

# ╔═╡ 7d293dcb-613d-4886-be67-71e15153c782
begin
	#parent = Crystal("NiPyC2_experiment.cif", convert_to_p1=false)
	parent = Crystal("test.cif", convert_to_p1=false)
	infer_bonds!(parent, true)
end

# ╔═╡ 6d118845-a8a8-463d-857f-aa9e9f49d795
view_structure(parent)

# ╔═╡ 5920a983-d74a-485c-bd4b-dddf95ed2b53


# ╔═╡ ee4a1efd-e90e-4298-8d57-56fec76c98d2
parent

# ╔═╡ 64df9de5-9bd0-40ec-8333-746cf26c09c8
query = moiety("3-H!-4-pyridyl.xyz")

# ╔═╡ b9731862-3472-4da1-a54c-88bd92facfc0
search = query in parent

# ╔═╡ 9c310fab-4c14-4a95-96ee-b3ed886c4584
temp = isomorphic_substructures(search)

# ╔═╡ cbdc90d4-fef1-4a48-a1a0-6c63c5711dbd
view_structure(temp)

# ╔═╡ 65fd8d82-cf66-4066-8d6f-c9b075f1e61e
view_structure(query)

# ╔═╡ 3debdd81-afdb-4cef-8c6d-3a0b80820965
replacement = moiety("3-F-4-pyridyl.xyz")

# ╔═╡ 34103ec4-8f0c-4e03-9c56-0c3e761d7692
view_structure(replacement)

# ╔═╡ 4f781ee2-dbf7-4676-ae41-82db98c20f7b
query in replacement

# ╔═╡ 8c980fa4-c5a8-4468-8a5f-c7902a38276c
primitive_child = replace(parent, query => replacement)

# ╔═╡ 1e4cce5a-9faa-409f-80aa-499e622fcdc2
write_cif(primitive_child, joinpath(rc[:paths][:crystals], "primitive_child.cif"))

# ╔═╡ 19217be5-c33b-4729-8259-d858a889f44b
write_xyz(primitive_child, joinpath(rc[:paths][:crystals], "primitive_child.xyz"))

# ╔═╡ 7d2d4404-3b29-4584-8e9a-69eb20cda3ef
replicated_child = Crystal("primitive_child.cif")

# ╔═╡ 8a6a4bef-5a29-4a53-998f-78748f4a0772
child = apply_symmetry_operations(primitive_child)

# ╔═╡ eb74a2bd-7e4f-48e6-adb9-4c6e2715df54
write_cif(child, joinpath(rc[:paths][:crystals], "replicated_child.cif"))

# ╔═╡ 0c79c15c-b807-4539-8389-97ad6f1bd102
npc = Crystal("NiPyC2_experiment.cif")

# ╔═╡ 85f1dc8b-75eb-42c8-9b90-833db362bdc8
infer_bonds!(npc, true)

# ╔═╡ 38c034fc-afba-446b-8878-39abfef39837
child2 = replace(npc, query => replacement)

# ╔═╡ 7493f873-7657-4535-a4d0-0d43035faef3
view_structure(child2)

# ╔═╡ d46787c6-59aa-4376-bd52-24395f90ca59
write_cif(child2, "child2.cif")

# ╔═╡ Cell order:
# ╠═9c36b2ed-28f1-492c-949d-4bf0e871152b
# ╠═42005e40-5450-11ec-3384-69ab8820daa2
# ╠═2b10d534-0c82-4c5f-8233-82ab7c519ce1
# ╠═7d293dcb-613d-4886-be67-71e15153c782
# ╠═b9731862-3472-4da1-a54c-88bd92facfc0
# ╠═9c310fab-4c14-4a95-96ee-b3ed886c4584
# ╠═cbdc90d4-fef1-4a48-a1a0-6c63c5711dbd
# ╠═6d118845-a8a8-463d-857f-aa9e9f49d795
# ╠═65fd8d82-cf66-4066-8d6f-c9b075f1e61e
# ╠═34103ec4-8f0c-4e03-9c56-0c3e761d7692
# ╠═5920a983-d74a-485c-bd4b-dddf95ed2b53
# ╠═4f781ee2-dbf7-4676-ae41-82db98c20f7b
# ╠═ee4a1efd-e90e-4298-8d57-56fec76c98d2
# ╠═64df9de5-9bd0-40ec-8333-746cf26c09c8
# ╠═3debdd81-afdb-4cef-8c6d-3a0b80820965
# ╠═8c980fa4-c5a8-4468-8a5f-c7902a38276c
# ╠═1e4cce5a-9faa-409f-80aa-499e622fcdc2
# ╠═19217be5-c33b-4729-8259-d858a889f44b
# ╠═7d2d4404-3b29-4584-8e9a-69eb20cda3ef
# ╠═8a6a4bef-5a29-4a53-998f-78748f4a0772
# ╠═eb74a2bd-7e4f-48e6-adb9-4c6e2715df54
# ╠═0c79c15c-b807-4539-8389-97ad6f1bd102
# ╠═85f1dc8b-75eb-42c8-9b90-833db362bdc8
# ╠═38c034fc-afba-446b-8878-39abfef39837
# ╠═7493f873-7657-4535-a4d0-0d43035faef3
# ╠═d46787c6-59aa-4376-bd52-24395f90ca59
