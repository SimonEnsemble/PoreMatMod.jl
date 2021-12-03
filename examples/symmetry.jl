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
	parent = Crystal("MOF-74.cif", convert_to_p1=false)
	infer_bonds!(parent, true)
	write_xyz(parent, "primitive_mof74.xyz")
end

# ╔═╡ 64df9de5-9bd0-40ec-8333-746cf26c09c8
query = moiety("primitive_fragment.xyz")

# ╔═╡ 3debdd81-afdb-4cef-8c6d-3a0b80820965
replacement = moiety("me_prim_frag.xyz")

# ╔═╡ 8c980fa4-c5a8-4468-8a5f-c7902a38276c
primitive_child = replace(parent, query => replacement)

# ╔═╡ 7ce201f3-b652-4e02-9c3b-b54bfc6f5bf8


# ╔═╡ 2beb47d7-f1ee-4ffe-adaf-9bccd2c3455f
view_structure(primitive_child)

# ╔═╡ b79321af-52f9-4c0f-ad52-0da265cd91df
primitive_child.symmetry

# ╔═╡ 19217be5-c33b-4729-8259-d858a889f44b
write_cif(primitive_child, joinpath(rc[:paths][:crystals], "primitive_child.cif"))

# ╔═╡ 7d2d4404-3b29-4584-8e9a-69eb20cda3ef
child = Crystal("primitive_child.cif")

# ╔═╡ Cell order:
# ╠═9c36b2ed-28f1-492c-949d-4bf0e871152b
# ╠═42005e40-5450-11ec-3384-69ab8820daa2
# ╠═2b10d534-0c82-4c5f-8233-82ab7c519ce1
# ╠═7d293dcb-613d-4886-be67-71e15153c782
# ╠═64df9de5-9bd0-40ec-8333-746cf26c09c8
# ╠═3debdd81-afdb-4cef-8c6d-3a0b80820965
# ╠═8c980fa4-c5a8-4468-8a5f-c7902a38276c
# ╠═7ce201f3-b652-4e02-9c3b-b54bfc6f5bf8
# ╠═2beb47d7-f1ee-4ffe-adaf-9bccd2c3455f
# ╠═b79321af-52f9-4c0f-ad52-0da265cd91df
# ╠═19217be5-c33b-4729-8259-d858a889f44b
# ╠═7d2d4404-3b29-4584-8e9a-69eb20cda3ef
