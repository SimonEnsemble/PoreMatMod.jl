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
	#parent = Crystal("NiPyC2_experiment.cif")
	parent = Crystal("test.cif")
	infer_bonds!(parent, true)
end

# ╔═╡ 6d118845-a8a8-463d-857f-aa9e9f49d795
view_structure(parent)

# ╔═╡ 64df9de5-9bd0-40ec-8333-746cf26c09c8
query = moiety("3-H!-4-pyridyl.xyz")

# ╔═╡ 65fd8d82-cf66-4066-8d6f-c9b075f1e61e
view_structure(query)

# ╔═╡ 3debdd81-afdb-4cef-8c6d-3a0b80820965
replacement = moiety("3-F-4-pyridyl.xyz")

# ╔═╡ 34103ec4-8f0c-4e03-9c56-0c3e761d7692
view_structure(replacement)

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

# ╔═╡ a1ab19c7-4aac-4c2f-a6b0-49492952acac
view_structure(child)

# ╔═╡ eb74a2bd-7e4f-48e6-adb9-4c6e2715df54
write_cif(child, joinpath(rc[:paths][:crystals], "replicated_child.cif"))

# ╔═╡ Cell order:
# ╠═9c36b2ed-28f1-492c-949d-4bf0e871152b
# ╠═42005e40-5450-11ec-3384-69ab8820daa2
# ╠═2b10d534-0c82-4c5f-8233-82ab7c519ce1
# ╠═7d293dcb-613d-4886-be67-71e15153c782
# ╠═6d118845-a8a8-463d-857f-aa9e9f49d795
# ╠═65fd8d82-cf66-4066-8d6f-c9b075f1e61e
# ╠═34103ec4-8f0c-4e03-9c56-0c3e761d7692
# ╠═64df9de5-9bd0-40ec-8333-746cf26c09c8
# ╠═3debdd81-afdb-4cef-8c6d-3a0b80820965
# ╠═8c980fa4-c5a8-4468-8a5f-c7902a38276c
# ╠═1e4cce5a-9faa-409f-80aa-499e622fcdc2
# ╠═19217be5-c33b-4729-8259-d858a889f44b
# ╠═7d2d4404-3b29-4584-8e9a-69eb20cda3ef
# ╠═8a6a4bef-5a29-4a53-998f-78748f4a0772
# ╠═a1ab19c7-4aac-4c2f-a6b0-49492952acac
# ╠═eb74a2bd-7e4f-48e6-adb9-4c6e2715df54
