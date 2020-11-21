### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ 2ae318e0-2c3f-11eb-2ceb-a1efd0374edd
# needed for development; remove when packages are registered
begin
	push!(LOAD_PATH, joinpath(homedir(), ".julia/dev/MOFfun.jl/src"))
	push!(LOAD_PATH, joinpath(homedir(), ".julia/dev/Xtals.jl/src"))
end;

# ╔═╡ e64622d0-2c3f-11eb-1872-21aab0c9aef8
using MOFun

# ╔═╡ d7b64d70-2c40-11eb-05c4-3d88c0c676cc
begin
	using PlutoUI, Bio3DView
	function viz_xtal(xtal::Crystal)
		write_xyz(xtal, joinpath(homedir(), ".mofun/temp/xyz.xyz"))
		write_vtk(xtal.box, joinpath(homedir(), ".mofun/temp/vtk.vtk"))
		viewfile(".mofun/temp/xyz.xyz", style=Style("stick"),
			vtkcell=".mofun/temp/vtk.vtk", axes=Axes(5, 0.25))
	end
end

# ╔═╡ a459dfc0-2c3e-11eb-1105-4d5377338a17
md"""
# MOFun.jl Tutorial Notebook

All `MOFun` uses need to start by [installing MOFun](https://github.com/SimonEnsemble/MOFun.jl) and loading the package into the namespace:
"""

# ╔═╡ ca1a6a70-2c40-11eb-0ca8-756b4e2317ca
md"""
For this notebook, we'll also need:
"""

# ╔═╡ b97f1830-2c42-11eb-0f5b-05751fa0b504
md"""
Tell `MOFun` where your data live with `set_path_to_data`:
"""

# ╔═╡ 14f5112e-2c46-11eb-2272-af228e1d0d20
set_path_to_data(joinpath(homedir(), ".mofun/data"), relpaths=true, print=true)

# ╔═╡ cbd4b2b0-2c42-11eb-3022-d57a58d90fc4
with_terminal() do 
	set_path_to_data(joinpath(homedir(), ".mofun/data"), relpaths=true, print=true)
end

# ╔═╡ 2089e850-2c40-11eb-0fc1-55038e192555
md"""
## Example: MOF Functionalization

Let's add the acetylamido group at the 2-position of the *p*-phenylene linker in IRMOF-1.

[illustration](illustration)
"""

# ╔═╡ 91a47530-2c42-11eb-16af-11e094126e3c
md"""
Load in the data files:
"""

# ╔═╡ 96753720-2c42-11eb-36f0-dd0ba2ad9bb5
begin
	IRMOF1 = Crystal("IRMOF-1.cif", infer_bonds=true)
	o_pphen = moiety("2-!-p-phenylene")
	a_pphen = moiety("2-acetylamido-p-phenylene")
end;

# ╔═╡ 968b5730-2c42-11eb-21f6-d16948182ab4
md"""
Run the search:
"""

# ╔═╡ 96a017b0-2c42-11eb-05a0-9be01e73d3e1
IRMOF1_pphen = substructure_search(o_pphen, IRMOF1);

# ╔═╡ 96b26730-2c42-11eb-1f98-5b337ad3b1a5
md"""
`Search` object:
"""

# ╔═╡ 96c83920-2c42-11eb-1b7a-c3aa70d41efa
IRMOF1_pphen.query.parent

# ╔═╡ e093af40-2c46-11eb-325b-73c18ab05b37
IRMOF1_pphen.query.s_moty

# ╔═╡ 96d95020-2c42-11eb-224a-c3995c16923e
IRMOF1_pphen.results

# ╔═╡ 96ee37b0-2c42-11eb-18c7-553effd43fff
md"""
Examine the results:
"""

# ╔═╡ 970198a0-2c42-11eb-053d-e9d0e12d97e7
nb_isomorphisms(IRMOF1_pphen) # how many total isomorphisms

# ╔═╡ 971520a0-2c42-11eb-38ef-d7dfe8951625
nb_locations(IRMOF1_pphen) # how many unique locations

# ╔═╡ 9728cfb0-2c42-11eb-299f-592bdff22fd5
nb_configs_at_loc(IRMOF1_pphen) # how many ways to orient a match at each location

# ╔═╡ 973c30a0-2c42-11eb-0508-cb0ae5fba58e
md"""
Make a random replacement at each location:
"""

# ╔═╡ 9752299e-2c42-11eb-23e7-1d2e3dbeb26b
begin
	xtal1 = find_replace(IRMOF1_pphen, a_pphen, rand_all=true)
	viz_xtal(xtal1)
end

# ╔═╡ 98355db0-2c42-11eb-1ae5-730f0d924589
md"""
Make a set number of random replacements:
"""

# ╔═╡ 987748ae-2c42-11eb-1cc0-b9df36b427c1
xtal2 = find_replace(IRMOF1_pphen, a_pphen, nb_loc=6)

# ╔═╡ 98b1b9a0-2c42-11eb-0202-17777eb4fa4a
md"""
Make a random replacement at selected locations:
"""

# ╔═╡ 99024aa0-2c42-11eb-17db-a12f2bb36d2d
xtal3 = find_replace(IRMOF1_pphen, a_pphen, loc=[2, 3, 5, 7])

# ╔═╡ 76c78e7e-2c44-11eb-0001-dd272ff02dd2
md"""
Make a specific series of replacements:
"""

# ╔═╡ 78618f70-2c44-11eb-2a2c-6581ac41c0de
xtal4 = find_replace(IRMOF1_pphen, a_pphen, loc=[4, 9, 16], ori=[1, 1, 1])

# ╔═╡ 78764ff0-2c44-11eb-21a8-a3e5e4944892


# ╔═╡ 7889b0e0-2c44-11eb-20b9-d9b6d5175837


# ╔═╡ 789fd0f0-2c44-11eb-3985-5963cf4ad4cf


# ╔═╡ 78b358f0-2c44-11eb-0dd0-15c887da4b07


# ╔═╡ 78c92ae0-2c44-11eb-0330-77e04414b29a


# ╔═╡ 78dcd9f0-2c44-11eb-10a1-25c29855c394


# ╔═╡ 78f061f0-2c44-11eb-35ae-675e4c6c407c


# ╔═╡ 79052270-2c44-11eb-1243-01f2779dca2e


# ╔═╡ 791b1b70-2c44-11eb-31e2-67a9cab69e6b


# ╔═╡ 792ea370-2c44-11eb-332b-a16df073c9c8


# ╔═╡ 79449c70-2c44-11eb-03f4-1b426a6f333f


# ╔═╡ 795a956e-2c44-11eb-2b07-3f922830a3d1


# ╔═╡ 796e1d70-2c44-11eb-2afc-0ddb4be3f2f1


# ╔═╡ 79841670-2c44-11eb-0bba-eb6cce97f686


# ╔═╡ 79979e70-2c44-11eb-0ec9-e3c6be14c0cf


# ╔═╡ 79ad9770-2c44-11eb-3044-ff0ea6f3ebad


# ╔═╡ 79c39072-2c44-11eb-0d47-fbec55ae49c5


# ╔═╡ Cell order:
# ╟─a459dfc0-2c3e-11eb-1105-4d5377338a17
# ╠═2ae318e0-2c3f-11eb-2ceb-a1efd0374edd
# ╠═e64622d0-2c3f-11eb-1872-21aab0c9aef8
# ╟─ca1a6a70-2c40-11eb-0ca8-756b4e2317ca
# ╠═d7b64d70-2c40-11eb-05c4-3d88c0c676cc
# ╟─b97f1830-2c42-11eb-0f5b-05751fa0b504
# ╠═14f5112e-2c46-11eb-2272-af228e1d0d20
# ╟─cbd4b2b0-2c42-11eb-3022-d57a58d90fc4
# ╟─2089e850-2c40-11eb-0fc1-55038e192555
# ╟─91a47530-2c42-11eb-16af-11e094126e3c
# ╠═96753720-2c42-11eb-36f0-dd0ba2ad9bb5
# ╟─968b5730-2c42-11eb-21f6-d16948182ab4
# ╠═96a017b0-2c42-11eb-05a0-9be01e73d3e1
# ╟─96b26730-2c42-11eb-1f98-5b337ad3b1a5
# ╠═96c83920-2c42-11eb-1b7a-c3aa70d41efa
# ╠═e093af40-2c46-11eb-325b-73c18ab05b37
# ╠═96d95020-2c42-11eb-224a-c3995c16923e
# ╟─96ee37b0-2c42-11eb-18c7-553effd43fff
# ╠═970198a0-2c42-11eb-053d-e9d0e12d97e7
# ╠═971520a0-2c42-11eb-38ef-d7dfe8951625
# ╠═9728cfb0-2c42-11eb-299f-592bdff22fd5
# ╟─973c30a0-2c42-11eb-0508-cb0ae5fba58e
# ╠═9752299e-2c42-11eb-23e7-1d2e3dbeb26b
# ╟─98355db0-2c42-11eb-1ae5-730f0d924589
# ╠═987748ae-2c42-11eb-1cc0-b9df36b427c1
# ╟─98b1b9a0-2c42-11eb-0202-17777eb4fa4a
# ╠═99024aa0-2c42-11eb-17db-a12f2bb36d2d
# ╟─76c78e7e-2c44-11eb-0001-dd272ff02dd2
# ╠═78618f70-2c44-11eb-2a2c-6581ac41c0de
# ╠═78764ff0-2c44-11eb-21a8-a3e5e4944892
# ╠═7889b0e0-2c44-11eb-20b9-d9b6d5175837
# ╠═789fd0f0-2c44-11eb-3985-5963cf4ad4cf
# ╠═78b358f0-2c44-11eb-0dd0-15c887da4b07
# ╠═78c92ae0-2c44-11eb-0330-77e04414b29a
# ╠═78dcd9f0-2c44-11eb-10a1-25c29855c394
# ╠═78f061f0-2c44-11eb-35ae-675e4c6c407c
# ╠═79052270-2c44-11eb-1243-01f2779dca2e
# ╠═791b1b70-2c44-11eb-31e2-67a9cab69e6b
# ╠═792ea370-2c44-11eb-332b-a16df073c9c8
# ╠═79449c70-2c44-11eb-03f4-1b426a6f333f
# ╠═795a956e-2c44-11eb-2b07-3f922830a3d1
# ╠═796e1d70-2c44-11eb-2afc-0ddb4be3f2f1
# ╠═79841670-2c44-11eb-0bba-eb6cce97f686
# ╠═79979e70-2c44-11eb-0ec9-e3c6be14c0cf
# ╠═79ad9770-2c44-11eb-3044-ff0ea6f3ebad
# ╠═79c39072-2c44-11eb-0d47-fbec55ae49c5
