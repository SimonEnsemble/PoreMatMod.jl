### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ f05016c0-3da9-11eb-3be7-01ad0b445dad
begin
	# needed for development. remove when packages registered
	push!(LOAD_PATH, joinpath(homedir(), ".julia/dev/MOFfun.jl/src"))
	push!(LOAD_PATH, joinpath(homedir(), ".julia/dev/Xtals.jl/src"))
end;

# ╔═╡ df039260-3daa-11eb-3434-0f02c2270267
using MOFun, Bio3DView, Xtals

# ╔═╡ 0ef8e7a0-3daa-11eb-2940-5fd52143a505
md"""
# Removing Adsorbed Solvent
Experimental structures of as-synthesized MOFs commonly contain solvent molecules which are stably adsorbed in the pores.  Removal of these solvent molecules is known as "activation."  Here is an example of how to activate a MOF structure *in silico*.
"""

# ╔═╡ ddd79580-3daa-11eb-1d0c-dd15eb627f6d
md"""
## Preparing Data
"""

# ╔═╡ a029a2f0-3daa-11eb-24dd-6dced09f460f
md"""
Obtain the structure of a [MOF with adsorbed solvent]() and the [solvent molecule]().
"""

# ╔═╡ a0d85980-3daa-11eb-31c9-07196ed505e9
md"""
## Loading Data
"""

# ╔═╡ ecfa5c10-3da9-11eb-3981-6f5ac224e7da
xtal = Crystal("MOF_w_solvent.cif")

# ╔═╡ 7a66e370-3daa-11eb-30a7-97dcbab917d3
infer_bonds!(xtal, true)

# ╔═╡ 7a677fb0-3daa-11eb-3618-21d74ba799ed
s_moty = moiety("solvent")

# ╔═╡ 7a697b80-3daa-11eb-0e53-719420d92111
r_moty = moiety(nothing)

# ╔═╡ afa25c40-3daa-11eb-01e9-b15600030792
md"""
## Find and Replace
"""

# ╔═╡ 7a6a17c0-3daa-11eb-2155-bbc7fdc88a08
activated_xtal = find_replace(s_moty ∈ xtal, r_moty, rand_all=true)

# ╔═╡ 007d4bd0-3dd2-11eb-382b-d5529fe3bdf4
begin
	no_pb = deepcopy(activated_xtal)
	Xtals.drop_cross_pb_bonds!(no_pb)
	Xtals.write_mol2(no_pb, filename="view.mol2")
	viewfile("view.mol2", "mol2")
end

# ╔═╡ Cell order:
# ╠═f05016c0-3da9-11eb-3be7-01ad0b445dad
# ╟─0ef8e7a0-3daa-11eb-2940-5fd52143a505
# ╠═df039260-3daa-11eb-3434-0f02c2270267
# ╟─ddd79580-3daa-11eb-1d0c-dd15eb627f6d
# ╟─a029a2f0-3daa-11eb-24dd-6dced09f460f
# ╟─a0d85980-3daa-11eb-31c9-07196ed505e9
# ╠═ecfa5c10-3da9-11eb-3981-6f5ac224e7da
# ╠═7a66e370-3daa-11eb-30a7-97dcbab917d3
# ╠═7a677fb0-3daa-11eb-3618-21d74ba799ed
# ╠═7a697b80-3daa-11eb-0e53-719420d92111
# ╟─afa25c40-3daa-11eb-01e9-b15600030792
# ╠═7a6a17c0-3daa-11eb-2155-bbc7fdc88a08
# ╠═007d4bd0-3dd2-11eb-382b-d5529fe3bdf4
