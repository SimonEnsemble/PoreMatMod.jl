### A Pluto.jl notebook ###
# v0.12.14

using Markdown
using InteractiveUtils

# ╔═╡ 0b31daf0-32dd-11eb-168a-f92297478213
begin
	# needed for development. remove when package registered
	push!(LOAD_PATH, joinpath(homedir(), ".julia/dev/MOFfun.jl/src"))
	push!(LOAD_PATH, joinpath(homedir(), ".julia/dev/Xtals.jl/src"))
	;
end

# ╔═╡ 29583ff0-32de-11eb-12d6-bdae9ce336e6
using MOFun, Bio3DView

# ╔═╡ 98622730-32de-11eb-30cf-579b111fb718
md"""
# Replacing BDC in IRMOF-1 with an Acetylamido Derivative

This is the main example used throughout the documentation for `MOFun`.
"""

# ╔═╡ edb278d0-3332-11eb-20d0-cff6c491fa14
md"""
## Prepare Data
"""

# ╔═╡ f3efea22-3332-11eb-3a7d-9f09b1f557a3


# ╔═╡ 21c46660-32de-11eb-2a7c-c9ad4994a40e
md"""
## Load Data
"""

# ╔═╡ e5a00492-32de-11eb-2515-d5f42dce050b
begin
	xtal = Crystal("IRMOF-1.cif", infer_bonds=:cordero, periodic_boundaries=true)
	s_moty = moiety("2-!-p-phenylene")
	r_moty = moiety("2-acetylamido-p-phenylene")
end

# ╔═╡ 3d1404b0-32df-11eb-1763-993377fbca52
md"""
## Perform the search
"""

# ╔═╡ 44fd9f60-32df-11eb-0bf6-0d97cf428224
search = substructure_search(s_moty, xtal)

# ╔═╡ 4983006e-32df-11eb-187d-9f7165eb6e4e
md"""
## Make the Replacement
"""

# ╔═╡ 56ec13a0-32df-11eb-177f-1b2aa24b5880
new_xtal = find_replace(search, r_moty, rand_all=true)

# ╔═╡ ba8594b0-3332-11eb-2a78-451d1bc1ad3a


# ╔═╡ 5c079a30-32df-11eb-0214-e57bb41c1759


# ╔═╡ db785590-3332-11eb-0ed5-c9cb11abcff1


# ╔═╡ Cell order:
# ╠═0b31daf0-32dd-11eb-168a-f92297478213
# ╟─98622730-32de-11eb-30cf-579b111fb718
# ╠═29583ff0-32de-11eb-12d6-bdae9ce336e6
# ╠═edb278d0-3332-11eb-20d0-cff6c491fa14
# ╠═f3efea22-3332-11eb-3a7d-9f09b1f557a3
# ╠═21c46660-32de-11eb-2a7c-c9ad4994a40e
# ╠═e5a00492-32de-11eb-2515-d5f42dce050b
# ╠═3d1404b0-32df-11eb-1763-993377fbca52
# ╠═44fd9f60-32df-11eb-0bf6-0d97cf428224
# ╠═4983006e-32df-11eb-187d-9f7165eb6e4e
# ╠═56ec13a0-32df-11eb-177f-1b2aa24b5880
# ╠═ba8594b0-3332-11eb-2a78-451d1bc1ad3a
# ╠═5c079a30-32df-11eb-0214-e57bb41c1759
# ╠═db785590-3332-11eb-0ed5-c9cb11abcff1
