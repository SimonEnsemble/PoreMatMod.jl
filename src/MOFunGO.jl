### A Pluto.jl notebook ###
# v0.12.14

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
begin
    push!(LOAD_PATH, joinpath(homedir(), ".julia/dev/MOFfun.jl/src"))
    push!(LOAD_PATH, joinpath(homedir(), ".julia/dev/Xtals.jl/src"))
    using MOFun, PlutoUI, Bio3DView
    HOME = joinpath(homedir(), ".mofungo")
    set_path_to_data(joinpath(HOME, "data"), relpaths=true)
    dirs = ["", "temp", "data"]
    for dir in dirs
        dir = joinpath(HOME, dir)
        if !isdir(dir)
            mkdir(dir)
        end
    end
    cd(HOME)
    @bind load_inputs Button("Reset")
end

# ╔═╡ 6c1969e0-02f5-11eb-3fa2-09931a63b1ac
md"""
# [💠 MOFunGO 🍌](https://en.wikipedia.org/wiki/Mofongo)

This notebook interactively substitutes moieties within a `Crystal` using a modified implementation of Ullmann's algorithm to perform substructure searches and applying singular value decomposition to align fragments of the generated materials. Read the docs [here](https://github.com/SimonEnsemble/MOFunGO.jl).

See the original publication on MOFun.jl here: [(article)](http://localhost:1234) [(GitHub)](https://github.com/SimonEnsemble/MOFun.jl)


"""

# ╔═╡ 50269ffe-02ef-11eb-0614-f11975d991fe
begin load_inputs
    # input fields: s_moty, r_moty, xtal
    md"""
    ##### Input Files

    Parent Crystal $(@bind parent_crystal FilePicker())

    Search Moiety $(@bind search_moiety FilePicker())

    Replace Moiety $(@bind replace_moiety FilePicker())
    """
end

# ╔═╡ 33b1fb50-0f73-11eb-2ab2-9d2cb6c5a533
# write file input strings to files in temp directory
begin
    # dict for tracking load status of inputs
    isloaded = Dict([:r_moty => false, :s_moty => false, :parent => false])
    # r_moty loader
    if replace_moiety["data"] != UInt8[]
        write("$HOME/temp/r_moty.xyz", replace_moiety["data"])
        r_moty = moiety("temp/r_moty")
        isloaded[:r_moty] = true
    end
    # s_moty loader
    if search_moiety["data"] != UInt8[]
        write("$HOME/temp/s_moty.xyz", search_moiety["data"])
        s_moty = moiety("temp/s_moty")
        isloaded[:s_moty] = true
    end
    # xtal loader
    if parent_crystal["data"] != UInt8[]
        write("$HOME/temp/parent.cif", parent_crystal["data"])
        xtal = Crystal("$HOME/temp/parent.cif", check_overlap=false)
        strip_numbers_from_atom_labels!(xtal)
        infer_bonds!(xtal, true)
        isloaded[:parent] = true
    end
    # run search and display terminal message
    if isloaded[:s_moty] && isloaded[:parent]
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
            output = md"Number of locations $(@bind nb_loc Slider(1:nb_locations(search)))"
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
        elseif replace_mode == "random replacement at n random locations" && nb_loc > 0
            new_xtal = find_replace(search, r_moty, nb_loc=nb_loc)
        elseif replace_mode == "random replacement at specific locations" && loc ≠ []
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

# ╔═╡ 5918f770-103d-11eb-0537-81036bd3e675
begin
    if new_xtal_flag
        write_cif(new_xtal, "$HOME/temp/cif.cif")
        write_xyz(new_xtal, "$HOME/temp/xyz.xyz")
        write_vtk(new_xtal.box, "$HOME/temp/box.vtk")
        write_bond_information(new_xtal, "$HOME/temp/bonds.vtk")
        viewfile("temp/xyz.xyz", "xyz", vtkcell="temp/box.vtk", axes=Axes(4, 0.25))
    end
end

# ╔═╡ 31832e30-1054-11eb-24ed-219fd3e236a1
if new_xtal_flag
    download_cif = DownloadButton(read("temp/cif.cif"), "crystal.cif")
    download_box = DownloadButton(read("temp/box.vtk"), "unit_cell.vtk")
    download_xyz = DownloadButton(read("temp/xyz.xyz"), "atoms.xyz")
    download_bonds = DownloadButton(read("temp/bonds.vtk"), "bonds.vtk")
    md"""
    ### Output Files
    $download_cif
    $download_box
    $download_xyz
    $download_bonds
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
# ╟─5918f770-103d-11eb-0537-81036bd3e675
# ╟─31832e30-1054-11eb-24ed-219fd3e236a1
# ╟─5dc43a20-10b8-11eb-26dc-7fb98e9aeb1a
# ╟─90696d20-10b7-11eb-20b5-6174faeaf613
