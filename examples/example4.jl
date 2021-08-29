### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 143fa303-2ce1-471e-ab1f-09a77b88eb75
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ 37939a7a-0651-11ec-11c1-6b5ef0a19ec2
using PoreMatMod, Bio3DView

# ╔═╡ 8d523993-6e85-443a-9949-12030552b457
md"""
## PoreMatMod.jl Example 4
"""

# ╔═╡ 21bf7b62-1aef-45c6-9da4-db8d4d69604c
function view_structure(xtal::Crystal)
	write_vtk(xtal.box, "unit_cell.vtk")
	no_pb = deepcopy(xtal)
	drop_cross_pb_bonds!(no_pb)
	write_mol2(no_pb, filename="view.mol2")
	viewfile("view.mol2", "mol2", vtkcell="unit_cell.vtk", axes=Axes(4, 0.25))
end;

# ╔═╡ 5b71d14a-be80-4ac3-8983-62571d0d4e7d
md"""
### Generate Missing-Linker Defect

Example: create a new channel in UiO-66 via missing-linker defects and formate ion capping.
"""

# ╔═╡ 74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
begin
	parent = Crystal("UiO-66.cif")
	infer_bonds!(parent, true)
	query = moiety("BDC.xyz")
	replacement = moiety("formate_caps.xyz")
	with_defect = replace(parent, query => replacement; loc=[1, 3, 9])
	view_structure(with_defect)
end

# ╔═╡ Cell order:
# ╟─143fa303-2ce1-471e-ab1f-09a77b88eb75
# ╟─8d523993-6e85-443a-9949-12030552b457
# ╠═37939a7a-0651-11ec-11c1-6b5ef0a19ec2
# ╠═21bf7b62-1aef-45c6-9da4-db8d4d69604c
# ╟─5b71d14a-be80-4ac3-8983-62571d0d4e7d
# ╠═74aa19d2-b1a4-4333-9ff9-e6ea74e7d989
