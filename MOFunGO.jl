### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 50269ffe-02ef-11eb-0614-f11975d991fe
# libraries and stuff
begin
	# Pluto doesn't inherit environment variables, it seems. 😠
	push!(LOAD_PATH, joinpath("C:\\Users\\eahen\\.julia\\dev\\MOFfun.jl\\src"))
	using PorousMaterials, MOFun, Logging, PlutoUI, Bio3DView
	# Pluto also doesn't like to run from ~/.julia/ 😖
	@eval MOFun PATH_TO_MOIETIES="C:\\Users\\eahen\\.julia\\dev\\MOFfun.jl\\data\\moieties"
	# Main console gets output and it isn't color-coded 😢
	global_logger(ConsoleLogger(stdout, Logging.Info))
end;

# ╔═╡ 13459850-03c9-11eb-06dc-d91bd1826b28
# This is from https://gist.github.com/GiggleLiu/aff2af66a896cf8a05310b8ba66f540f#file-plutouitips-jl :
html"""
Expand/Contract Code <button id="showhide">↔</button>
<style>
	body.hide_all_inputs pluto-input {
		display: none;
	}
	body.hide_all_inputs pluto-shoulder {
		display: none;
	}
	body.hide_all_inputs pluto-trafficlight {
		display: none;
	}
	body.hide_all_inputs pluto-runarea {
		display: none;
	}
	body.hide_all_inputs .add_cell {
		display: none;
	}
	body.hide_all_inputs pluto-cell {
		min-height: 0;
		margin-top: 10px;
	}
</style>
<script>
	const button = this.querySelector("#showhide");
	button.onclick = () => {
		document.body.classList.toggle("hide_all_inputs")
	}
</script>
"""

# ╔═╡ 6c1969e0-02f5-11eb-3fa2-09931a63b1ac
# title/info
# Made with some inspiration from [https://gist.github.com/GiggleLiu](https://gist.github.com/GiggleLiu) 
md"""
# [💠 MOFunGO 🍌](https://en.wikipedia.org/wiki/Mofongo)

This notebook interactively substitutes moieties within a `Crystal` using a modified implementation of Ullmann's algorithm to perform substructure searches and applying the orthogonal Procrustes problem to align fragments of the generated materials.
"""

# ╔═╡ 77d70610-02ec-11eb-0004-fff39ffbb197
# input fields: s_moty, r_moty, xtal
md"""
### Find/Replace Options

Parent Crystal $(@bind PARENT TextField(default="IRMOF-1.cif"))

Search Moiety $(@bind SEARCH_MOIETY TextField(default="find-replace/2-!-p-phenylene"))

Replace Moiety $(@bind REPLACE_MOIETY TextField(default="2-acetylamido-p-phenylene"))
"""

# ╔═╡ 47ca486e-03c2-11eb-3b76-6337fc07c447
# vis refresh button
md"""
### Substitute and Visualize $(@bind UPDATE Button("GO"))
"""

# ╔═╡ 8768f630-03bc-11eb-0be1-338983860b0d
# refreshes visualization on button press
let UPDATE
	# not sure why this isn't rendering the unit cell, but it also doesn't in the tutorial notebook...
	try
		viewfile("MOFunGO_temp.xyz", "xyz", vtkcell="MOFunGO_temp.vtk")
	catch # suppress errors when no visuals have processed yet
	end
end

# ╔═╡ 41448e10-03c3-11eb-2a21-2102e629ffca
# file output
md"""
#### Save Result
$(@bind SAVE_PATH TextField(default="new_xtal"))
.cif $(@bind SAVE_CIF CheckBox())
.vtk $(@bind SAVE_VTK CheckBox())
$(@bind SAVE_RESULT Button("Save"))
"""

# ╔═╡ c4865f10-02ec-11eb-3cb1-8fa065c4b09f
begin
	xtal = Crystal(PARENT)
	strip_numbers_from_atom_labels!(xtal)
	infer_bonds!(xtal, true)
	s_moty = moiety(SEARCH_MOIETY)
	r_moty = moiety(REPLACE_MOIETY)
end;

# ╔═╡ 685ff86e-0302-11eb-06fe-11a418009951
begin
	search = s_moty ∈ xtal
end;

# ╔═╡ 5fff5d30-02fb-11eb-011c-e9af7248c73e
begin
	locations = ["$(loc)" for loc in 1:search.num_locations]
end;

# ╔═╡ fea0b040-03c0-11eb-2fc4-cb02d90a5be1
# location selector
md"""
Location in Parent $(@bind LOCATION Select(locations))
"""

# ╔═╡ 76bf2b6e-0302-11eb-2331-8dac63ed0b94
begin
	orientations = ["$(res.configuration.orientation)" for res in search.results if "$(res.configuration.location)" == LOCATION]
end;

# ╔═╡ f705aac2-03c0-11eb-1185-eb29469be398
# orientation selector
md"""
Orientation at Location $(@bind ORIENTATION Select(orientations))
"""

# ╔═╡ f723c59e-03b6-11eb-3e5d-15aa0edd28c3
begin
	new_xtal = find_replace(s_moty, r_moty, xtal, config=Configuration(parse(Int, LOCATION), parse(Int, ORIENTATION)))
end;

# ╔═╡ d8c9faf0-03bd-11eb-15b5-7b74057459c5
# terminal output
begin
	with_terminal() do
		@info new_xtal
	end
end

# ╔═╡ 57d03120-03c3-11eb-0b53-67c02471b008
let SAVE_RESULT
	if SAVE_CIF
		write_cif(new_xtal, SAVE_PATH*".cif")
		
	end
	if SAVE_VTK
		write_vtk(new_xtal.box, SAVE_PATH*".vtk")
	end
end

# ╔═╡ 06eaee10-03ac-11eb-0743-8f9c7bdfea5a
begin
	write_xyz(new_xtal, "MOFunGO_temp.xyz")
	write_vtk(new_xtal.box, "MOFunGO_temp.vtk")
end
# visual doesn't update w/o button, but when button exists, no need to click it! (???)

# ╔═╡ 3b0134f0-02f5-11eb-0237-15a21c849f2a
# This is from https://gist.github.com/GiggleLiu/aff2af66a896cf8a05310b8ba66f540f#file-plutouitips-jl :
html"""
Expand/Contract Code <button id="showhide">↔</button>
<style>
	body.hide_all_inputs pluto-input {
		display: none;
	}
	body.hide_all_inputs pluto-shoulder {
		display: none;
	}
	body.hide_all_inputs pluto-trafficlight {
		display: none;
	}
	body.hide_all_inputs pluto-runarea {
		display: none;
	}
	body.hide_all_inputs .add_cell {
		display: none;
	}
	body.hide_all_inputs pluto-cell {
		min-height: 0;
		margin-top: 10px;
	}
</style>
<script>
	const button = this.querySelector("#showhide");
	button.onclick = () => {
		document.body.classList.toggle("hide_all_inputs")
	}
</script>
"""

# ╔═╡ Cell order:
# ╠═13459850-03c9-11eb-06dc-d91bd1826b28
# ╠═6c1969e0-02f5-11eb-3fa2-09931a63b1ac
# ╠═50269ffe-02ef-11eb-0614-f11975d991fe
# ╠═77d70610-02ec-11eb-0004-fff39ffbb197
# ╠═fea0b040-03c0-11eb-2fc4-cb02d90a5be1
# ╠═f705aac2-03c0-11eb-1185-eb29469be398
# ╠═47ca486e-03c2-11eb-3b76-6337fc07c447
# ╠═8768f630-03bc-11eb-0be1-338983860b0d
# ╠═d8c9faf0-03bd-11eb-15b5-7b74057459c5
# ╠═41448e10-03c3-11eb-2a21-2102e629ffca
# ╠═57d03120-03c3-11eb-0b53-67c02471b008
# ╠═c4865f10-02ec-11eb-3cb1-8fa065c4b09f
# ╠═685ff86e-0302-11eb-06fe-11a418009951
# ╠═5fff5d30-02fb-11eb-011c-e9af7248c73e
# ╠═76bf2b6e-0302-11eb-2331-8dac63ed0b94
# ╠═06eaee10-03ac-11eb-0743-8f9c7bdfea5a
# ╠═f723c59e-03b6-11eb-3e5d-15aa0edd28c3
# ╠═3b0134f0-02f5-11eb-0237-15a21c849f2a
