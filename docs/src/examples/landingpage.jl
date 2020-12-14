### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 6e9609e0-32dd-11eb-2ae9-5b71bd231df6
begin
	# needed for development. remove when packages registered
	push!(LOAD_PATH, joinpath(homedir(), ".julia/dev/MOFfun.jl/src"))
	push!(LOAD_PATH, joinpath(homedir(), ".julia/dev/Xtals.jl/src"))
end;

# ╔═╡ d22fba50-3daa-11eb-2db8-435563bce680
using MOFun, Bio3DView

# ╔═╡ da992d02-3da2-11eb-2b1e-73394f0feaf0
md"""
# Repairing Disorder, Activating the MOF, and Derivatizing
Here is a demonstration of a multi-step work flow.
"""

# ╔═╡ 367d0a60-3da3-11eb-2831-7342dfff73a2
md"""
## Preparing Data
"""

# ╔═╡ 39cbe740-3da3-11eb-3b5a-b33fc04239d6
md"""
For the starting point, obtain the [EMEHUB]() record from [CSD](https://www.ccdc.cam.ac.uk/solutions/csd-core/components/csd/).  Repairing the disorder requires extracting the [disordered PyC2 ligand]() and providing the [discrete replacement]().  Removal of the acetylene guest molecules requires [their structure]().  Finally, methylation of the PyC2 ligand requires the structure of the [derivatized ligand]().
"""

# ╔═╡ 3a6eb6f0-3da3-11eb-1be2-576825f3d690
md"""
## Loading Data
"""

# ╔═╡ 40fa965e-3da3-11eb-3ca9-bd2ef48129ce
xtal = Crystal("EMEHUB.cif", infer_bonds=:voronoi, periodic_boundaries=true)

# ╔═╡ 410dd040-3da3-11eb-3a49-17e008b3cc07
md"""
## Perform the Search
"""

# ╔═╡ 42fdb960-3da3-11eb-387c-293f474e02a0


# ╔═╡ 43102ff2-3da3-11eb-0a00-b90ab8d5e64f
md"""
## Make the Replacement
"""

# ╔═╡ f6f702d0-3daa-11eb-1df3-3191fb658738
begin
	
	repaired = (moiety("disordered_bipy!") => moiety("discrete")) ∈ xtal
	active = (moiety("acetylene") => moiety(nothing)) ∈ repaired
	novel = (moiety("2-H!-PyC2") => moiety("2-Me-PyC2")) ∈ active
end

# ╔═╡ Cell order:
# ╠═6e9609e0-32dd-11eb-2ae9-5b71bd231df6
# ╟─da992d02-3da2-11eb-2b1e-73394f0feaf0
# ╠═d22fba50-3daa-11eb-2db8-435563bce680
# ╟─367d0a60-3da3-11eb-2831-7342dfff73a2
# ╟─39cbe740-3da3-11eb-3b5a-b33fc04239d6
# ╟─3a6eb6f0-3da3-11eb-1be2-576825f3d690
# ╠═40fa965e-3da3-11eb-3ca9-bd2ef48129ce
# ╟─410dd040-3da3-11eb-3a49-17e008b3cc07
# ╠═42fdb960-3da3-11eb-387c-293f474e02a0
# ╟─43102ff2-3da3-11eb-0a00-b90ab8d5e64f
# ╠═f6f702d0-3daa-11eb-1df3-3191fb658738
