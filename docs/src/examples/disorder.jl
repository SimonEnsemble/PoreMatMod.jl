### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 5586f040-32dd-11eb-0f22-758ecc5823f1
begin
	# needed for development. remove when packages registered
	push!(LOAD_PATH, joinpath(homedir(), ".julia/dev/MOFfun.jl/src"))
	push!(LOAD_PATH, joinpath(homedir(), ".julia/dev/Xtals.jl/src"))
end;

# ╔═╡ cf1497a0-3daa-11eb-1a8e-1b397b937b01
using MOFun, Bio3DView, Xtals

# ╔═╡ adbf16f0-3da2-11eb-105f-5101165f4a1c
md"""
# Repairing Disorder
This is a common task when working with experimental data.  X-ray crystallography sometimes determines certain atoms to be statistically distributed among certain physical locations over all unit cells, and these atoms may appear together in the final spatial representation.  Here is how to replace a disordered moiety with its discrete representation.
"""

# ╔═╡ e34617b0-3da2-11eb-2870-e7ed50511ef5
md"""
## Preparing Data
"""

# ╔═╡ e7e3e2c0-3da2-11eb-0ccb-7bcaebee1bf1
md"""
Obtain the experimental structure of [ZmID](), extract its [disordered DABCO]() moiety, and prepare a [replacement DABCO]().
"""

# ╔═╡ e7f4d2b2-3da2-11eb-2b23-3b82aa2f99bf
md"""
## Loading Data
"""

# ╔═╡ f13d3a10-3da2-11eb-3aeb-8d20b4566755
begin
	xtal = Crystal("ZmID.cif", check_overlap=false)
	infer_bonds!(xtal, true)
	s_moty = moiety("disordered_dabco")
	r_moty = moiety("dabco")

end

# ╔═╡ f14b43d0-3da2-11eb-06a7-e579951c7ef1
md"""
## Perform the Search
"""

# ╔═╡ f701d1e0-3da2-11eb-12d3-eb05b282c507


# ╔═╡ f70ef140-3da2-11eb-0725-39af145ec116
md"""
## Make the Replacement
"""

# ╔═╡ f3fcfcb0-3daa-11eb-143e-d3744604602d
repaired_xtal = replace(s_moty ∈ xtal, r_moty, rand_all=true)

# ╔═╡ 1e77b440-3dd2-11eb-04cc-0daef8c13b8a
begin
	no_pb = deepcopy(repaired_xtal)
	Xtals.drop_cross_pb_bonds!(no_pb)
	Xtals.write_mol2(no_pb, filename="view.mol2")
	viewfile("view.mol2", "mol2")
end

# ╔═╡ Cell order:
# ╠═5586f040-32dd-11eb-0f22-758ecc5823f1
# ╟─adbf16f0-3da2-11eb-105f-5101165f4a1c
# ╠═cf1497a0-3daa-11eb-1a8e-1b397b937b01
# ╟─e34617b0-3da2-11eb-2870-e7ed50511ef5
# ╠═e7e3e2c0-3da2-11eb-0ccb-7bcaebee1bf1
# ╟─e7f4d2b2-3da2-11eb-2b23-3b82aa2f99bf
# ╠═f13d3a10-3da2-11eb-3aeb-8d20b4566755
# ╟─f14b43d0-3da2-11eb-06a7-e579951c7ef1
# ╠═f701d1e0-3da2-11eb-12d3-eb05b282c507
# ╟─f70ef140-3da2-11eb-0725-39af145ec116
# ╠═f3fcfcb0-3daa-11eb-143e-d3744604602d
# ╠═1e77b440-3dd2-11eb-04cc-0daef8c13b8a
