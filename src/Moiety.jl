module Moiety
export moiety, filter_R_group!, PATH_TO_MOIETIES

PATH_TO_MOIETIES = joinpath(pwd(), "data/moieties")

using PorousMaterials, DataFrames, LightGraphs


"""
Removes trailing underscores from the atoms of the input Crystal and returns their indices
"""
function filter_R_group!(xtal::Crystal)::Array{Int}
	@debug "Filtering R group"
	R = []
	for (idx, label) in enumerate(xtal.atoms.species) # loop over crystal atoms to find tags
		# if String representation of label Symbol ends in _, atom is in R
		tokens = split("$label", '_')
		if length(tokens) > 1 && tokens[2] == "" # other _ in symbol may be problematic.
			push!(R, idx)
			xtal.atoms.species[idx] = Symbol(tokens[1])
		end
	end
	@debug "Returning" R
	return R
end


"""
Tags atoms of specified indices with a trailing _ using df.index order
"""
function tag_R_group!(xtal::Crystal, arr::Array{Int}, df::DataFrame)
    # for each index in arr, find its row in df and change that label in xtal
    for r in arr # loop over R array
        i = getindex(df.index, r) # Map R label to new node order
        label = xtal.atoms.species[i] # Get label to edit
        xtal.atoms.species[i] = Symbol("$(label)_") # tag atom
    end
end


@doc raw"""
Generates a moiety (Crystal) from an .xyz file found in PATH_TO_MOIETIES
"""
function moiety(name::String)::Crystal
	@debug "Getting moiety: $name"
	# generate Crystal from moiety XYZ coords
    box = unit_cube()
    fx = Frac(read_xyz(joinpath(pwd(), "$PATH_TO_MOIETIES/$name.xyz")), box)
    charges = Charges{Frac}(0)
	moiety = Crystal(name, box, fx, charges)
	R_group = filter_R_group!(moiety) # collect R indices and un-tag atoms for bonding

	# sort by node degree (only needed for search moiety, but hurts nothing)
	infer_bonds!(moiety, false)
	df = DataFrame([[1:nv(moiety.bonds)...], degree(moiety.bonds)], [:index, :degree])
	sort!(df, :degree, rev=true)

	# ordered atoms w/o R group
	not_R = [i for i in 1:length(df.index) if ! (i in R_group)]
	order_wo_R = df.index[not_R]
	# append R-group to the end
	order = vcat(order_wo_R, R_group)
	@debug order moiety.atoms.species[order], moiety.atoms.coords[order]
	R_group = [length(order):-1:(length(order) - length(R_group) + 1)...]

	# rebuild Atoms
	atoms = Atoms(moiety.atoms.species[order], moiety.atoms.coords[order]) # atoms are sorted by degree and un-tagged

	# build moiety with ordered atoms
	moiety = Crystal(name, box, atoms, charges)
    infer_bonds!(moiety, false)
	if length(R_group) ≠ 0
		tag_R_group!(moiety, R_group, df) # replace tags
	end

	return moiety # nodes are sorted by bond order, and R group is tagged w/ _
end

end # module
