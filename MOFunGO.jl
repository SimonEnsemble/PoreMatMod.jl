### A Pluto.jl notebook ###
# v0.12.4

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

# ╔═╡ 90696d20-10b7-11eb-20b5-6174faeaf613
# gist.github.com/GiggleLiu/aff2af66a896cf8a05310b8ba66f540f#file-plutouitips-jl
# used to be able to toggle collapsing/expanding of all code in the notebook...
begin
	push!(LOAD_PATH, joinpath(homedir(), ".julia/dev/MOFfun.jl/src"))
	using PorousMaterials, MOFun, PlutoUI, Bio3DView
	@eval PorousMaterials PATH_TO_DATA=joinpath(homedir(), ".mofungo/data")
	@eval MOFun PATH_TO_MOIETIES=joinpath(PorousMaterials.PATH_TO_DATA, "moieties")
	HOME = joinpath(homedir(), ".mofungo")
	dirs = ["", "temp", "data", "data/moieties", "data/crystals"]
	for dir in dirs
		dir = joinpath(HOME, dir)
		if !isdir(dir)
			mkdir(dir)
		end
	end
	@bind load_inputs Button("Reset")
end

# ╔═╡ 6c1969e0-02f5-11eb-3fa2-09931a63b1ac
md"""
# [💠 MOFunGO 🍌](https://en.wikipedia.org/wiki/Mofongo)

This notebook interactively substitutes moieties within a `Crystal` using a modified implementation of Ullmann's algorithm to perform substructure searches and applying singular value decomposition to align fragments of the generated materials.

See the original publication on MOFun.jl here: [(article)](http://localhost:1234) [(GitHub)](http://localhost:1234)


"""

# ╔═╡ 50269ffe-02ef-11eb-0614-f11975d991fe
begin load_inputs
	# input fields: s_moty, r_moty, xtal
	md"""
	### Structure Inputs

	Parent Crystal $(@bind parent_crystal FilePicker()))

	Search Moiety $(@bind search_moiety FilePicker()))

	Replace Moiety $(@bind replace_moiety FilePicker()))
	"""
end

# ╔═╡ 33b1fb50-0f73-11eb-2ab2-9d2cb6c5a533
# write file input strings to files in temp directory
begin
	# dict for tracking load status of inputs
	isloaded = Dict([:r_moty => false, :s_moty => false, :parent => false])
	# r_moty loader
	if replace_moiety["data"] != UInt8[]
		write("$HOME/data/moieties/r_moty.xyz", replace_moiety["data"])
		r_moty = moiety("r_moty")
		isloaded[:r_moty] = true
	end
	# s_moty loader
	if search_moiety["data"] != UInt8[]
		write("$HOME/data/moieties/s_moty.xyz", search_moiety["data"])
		s_moty = moiety("s_moty")
		isloaded[:s_moty] = true
	end
	# xtal loader
	if parent_crystal["data"] != UInt8[]
		write("$HOME/data/parent.cif", parent_crystal["data"])
		xtal = Crystal("$HOME/temp/parent.cif")
		strip_numbers_from_atom_labels!(xtal)
		infer_bonds!(xtal, true)
		isloaded[:parent] = true
	end
	# run search and display terminal message
	if all(values(isloaded))
		search = s_moty ∈ xtal
		with_terminal() do
			@info "Search Results" isomorphisms=nb_isomorphisms(search) locations=nb_locations(search)
		end
	end
end

# ╔═╡ 415e9210-0f71-11eb-15c8-e7484b5be309
# choose replacement type
if all(values(isloaded))
	md"""
	### Find/Replace Options

	Mode $(@bind replace_mode Select(["", "random replacement at each location", "random replacement at n random locations", "random replacement at specific locations", "specific replacements"]))
	"""
end

# ╔═╡ 3997c4d0-0f75-11eb-2976-c161879b8d0c
# options populated w/ conditional logic based on mode selection
begin
	local output = nothing
	if all(values(isloaded))
		x = ["$(x)" for x in 1:nb_locations(search)]
		if replace_mode == "random replacement at each location"
			output = nothing
		elseif replace_mode == "random replacement at n random locations"
			output = md"Number of locations $(@bind nb_loc Slider(0:nb_locations(search)))"
		elseif replace_mode == "random replacement at specific locations"
			output = md"Locations $(@bind loc MultiSelect(x))"
		elseif replace_mode == "specific replacements"
			output = md"""
		Locations $(@bind loc MultiSelect(x))

		Orientations $(@bind ori TextField())
		"""
		else
			output = nothing
		end
		output
	end
end

# ╔═╡ 69edca20-0f94-11eb-13ba-334438ca2406
begin
	new_xtal_flag = false
	if all(values(isloaded))
		new_xtal_flag = true
		if replace_mode == "random replacement at each location"
			new_xtal = find_replace(search, r_moty, rand_all=true)
		elseif replace_mode == "random replacement at n random locations"
			new_xtal = find_replace(search, r_moty, nb_loc=nb_loc)
		elseif replace_mode == "random replacement at specific locations"
			new_xtal = find_replace(search, r_moty, loc=[parse(Int, x) for x in loc])
		elseif replace_mode == "specific replacements"
			if loc ≠ [] && ori ≠ "" && length(loc) == length(split(ori, ","))
				new_xtal = find_replace(search, r_moty, 
					loc=[parse(Int, x) for x in loc],
					ori=[parse(Int, x) for x in split(ori, ",")])
			else
				new_xtal_flag = false
			end
		else
			new_xtal_flag = false
		end
		if new_xtal_flag
			with_terminal() do
				if replace_mode == "random replacement at each location"
					@info replace_mode new_xtal
				elseif replace_mode == "random replacement at n random locations"
					@info replace_mode nb_loc new_xtal
				elseif replace_mode == "random replacement at specific locations"
					@info replace_mode loc new_xtal
				elseif replace_mode == "specific replacements"
					@info replace_mode loc ori new_xtal
				end
			end
		end
	end
end

# ╔═╡ 84a3724e-1059-11eb-3586-8dbef7b1233f
if new_xtal_flag
	md"""
	### Visualization
	
	(This is broken now, and I don't know why.)
	
	Enable $(@bind vis_czbx CheckBox())
	"""
end

# ╔═╡ 5918f770-103d-11eb-0537-81036bd3e675
if new_xtal_flag
	write_cif(new_xtal, "new_xtal.cif")
	write_vtk(new_xtal.box, "new_xtal")
	vis_czbx ? viewfile("new_xtal.cif", "cif", style=Style("stick"), vtkcell="new_xtal.vtk", axes=Axes(4, 0.25)) : nothing
end

# ╔═╡ 31832e30-1054-11eb-24ed-219fd3e236a1
if new_xtal_flag
	md"""
	### Output Files
	$(DownloadButton(read("new_xtal.cif"), "MOFunGO.cif")) 
	$(DownloadButton(read("new_xtal.vtk"), "MOFunGO.vtk"))
	"""
end

# ╔═╡ 5dc43a20-10b8-11eb-26dc-7fb98e9aeb1a
md"""
Adrian Henle, [Simon Ensemble](http://simonensemble.github.io), 2020

$(Resource("https://simonensemble.github.io/osu_logo.jpg", :width => 250))
"""

# ╔═╡ Cell order:
# ╟─6c1969e0-02f5-11eb-3fa2-09931a63b1ac
# ╟─50269ffe-02ef-11eb-0614-f11975d991fe
# ╟─33b1fb50-0f73-11eb-2ab2-9d2cb6c5a533
# ╟─415e9210-0f71-11eb-15c8-e7484b5be309
# ╟─3997c4d0-0f75-11eb-2976-c161879b8d0c
# ╟─69edca20-0f94-11eb-13ba-334438ca2406
# ╟─84a3724e-1059-11eb-3586-8dbef7b1233f
# ╟─5918f770-103d-11eb-0537-81036bd3e675
# ╟─31832e30-1054-11eb-24ed-219fd3e236a1
# ╟─90696d20-10b7-11eb-20b5-6174faeaf613
# ╟─5dc43a20-10b8-11eb-26dc-7fb98e9aeb1a
