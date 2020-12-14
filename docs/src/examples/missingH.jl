### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 64536400-32dd-11eb-3595-2d3665b19931
begin
	# needed for development. remove when packages registered
	push!(LOAD_PATH, joinpath(homedir(), ".julia/dev/MOFfun.jl/src"))
	push!(LOAD_PATH, joinpath(homedir(), ".julia/dev/Xtals.jl/src"))
end;

# ╔═╡ d4f40a6e-3daa-11eb-0574-5da8c9051ab5
using MOFun, Bio3DView, Xtals

# ╔═╡ 060b75ae-3da3-11eb-03c9-15d146120d6c
md"""
# Inserting Missing Hydrogens
Experimental structure data frequently omits the locations of hydrogen atoms, as they are difficult to solve precisely unless the experimental sample is of extremely good quality.  Here is how to addmissing  hydrogens back into a structure.
"""

# ╔═╡ 47c10dd0-3da3-11eb-09b7-d9c883f4a381
md"""
## Preparing Data
"""

# ╔═╡ 49e0e08e-3da3-11eb-274b-1bc81bde0bee
md"""
Obtain the structure of IRMOF-1 with [missing hydrogens](), extract the [tetra-dehydro-*p*-phenylene](), and prepare the properly hydrogenated [replacement moiety]().
"""

# ╔═╡ 49f63d50-3da3-11eb-088a-d1914af4dbe0
md"""
## Loading Data
"""

# ╔═╡ 4bcf6a20-3da3-11eb-04cc-71095cbc8e09
begin
	xtal = Crystal("IRMOF-1_noH.cif")
	s_moty = moiety("p-phenylene_noH")
	r_moty = moiety("p-phenylene")
end

# ╔═╡ 4be22ed0-3da3-11eb-1c60-4dd719c5f7f1
md"""
## Perform the Search
"""

# ╔═╡ 4d4408c0-3da3-11eb-2691-efa4cd9adc6d


# ╔═╡ 4d582d00-3da3-11eb-2cbd-c381e3183e10
md"""
## Make the Replacement
"""

# ╔═╡ 28be3c72-3dab-11eb-257f-cd33683a6fa3
repaired_xtal = find_replace(s_moty ∈ xtal, r_moty, rand_all=true)

# ╔═╡ 6f89cb20-3dd2-11eb-2de0-79aba5d6c6f5
begin
	no_pb = deepcopy(repaired_xtal)
	Xtals.drop_cross_pb_bonds!(no_pb)
	Xtals.write_mol2(no_pb, filename="view.mol2")
	viewfile("view.mol2", "mol2")
end

# ╔═╡ Cell order:
# ╠═64536400-32dd-11eb-3595-2d3665b19931
# ╟─060b75ae-3da3-11eb-03c9-15d146120d6c
# ╠═d4f40a6e-3daa-11eb-0574-5da8c9051ab5
# ╟─47c10dd0-3da3-11eb-09b7-d9c883f4a381
# ╠═49e0e08e-3da3-11eb-274b-1bc81bde0bee
# ╟─49f63d50-3da3-11eb-088a-d1914af4dbe0
# ╠═4bcf6a20-3da3-11eb-04cc-71095cbc8e09
# ╟─4be22ed0-3da3-11eb-1c60-4dd719c5f7f1
# ╠═4d4408c0-3da3-11eb-2691-efa4cd9adc6d
# ╟─4d582d00-3da3-11eb-2cbd-c381e3183e10
# ╠═28be3c72-3dab-11eb-257f-cd33683a6fa3
# ╠═6f89cb20-3dd2-11eb-2de0-79aba5d6c6f5
