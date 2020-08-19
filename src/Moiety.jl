module Moiety
export moiety, filter_R_group, PATH_TO_MOIETIES

PATH_TO_MOIETIES = joinpath(pwd(), "data/moieties")
R_GROUP_TAG = '!'

using PorousMaterials, DataFrames, LightGraphs


"""
Returns bonding rules including R-group-tagged atom copies
"""
function new_bonding_rules()::Array{BondingRule}
	bondingrules = PorousMaterials.default_bondingrules()
	push!(bondingrules, BondingRule(:C, :*, 0.4, 1.9))
	newrules = []
	# this loop is stupidly structured because the non-stupid version would hang, inexplicably, forever
	for i in 1:length(bondingrules)
	    if bondingrules[i].species_i != :*
	        push!(newrules, BondingRule(Symbol("$(bondingrules[i].species_i)!"), bondingrules[i].species_j, bondingrules[i].min_dist, bondingrules[i].max_dist))
	        push!(newrules, BondingRule(Symbol("$(bondingrules[i].species_j)!"), bondingrules[i].species_i, bondingrules[i].min_dist, bondingrules[i].max_dist))
	        push!(newrules, BondingRule(Symbol("$(bondingrules[i].species_i)!"), Symbol("$(bondingrules[i].species_j)!"), bondingrules[i].min_dist, bondingrules[i].max_dist))
	    end
	end
	return unique(vcat(bondingrules, newrules))
end


"""
Returns R group indices
"""
function filter_R_group(xtal::Crystal; remove = false)::Array{Int}
	@debug "Filtering R group" R_GROUP_TAG
	R = []
	for (idx, label) in enumerate(xtal.atoms.species) # loop over crystal atoms to find tags
		# if String representation of label Symbol ends in !, atom is in R
		tokens = split("$label", R_GROUP_TAG)
		if length(tokens) == 2 && tokens[2] == "" # other ! in symbol not tolerated.
			push!(R, idx)
			if remove
				xtal.atoms.species[idx] = Symbol(tokens[1])
			end
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
        xtal.atoms.species[i] = Symbol("$(label)%") # tag atom
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
	R_group_indices = filter_R_group!(moiety) # collect R indices and un-tag atoms for bonding

	# sort by node degree (only needed for search moiety, but hurts nothing)
	infer_bonds!(moiety, false)
	df = DataFrame([[1:nv(moiety.bonds)...], degree(moiety.bonds)], [:index, :degree])
	sort!(df, :degree, rev=true)
	## TODO change this to sortperm() for conciseness

	# ordered atoms w/o R group
	not_R = [i for i in 1:length(df.index) if ! (i in R_group_indices)]
	order_wo_R = df.index[not_R]
	# append R-group to the end
	order = vcat(order_wo_R, R_group_indices)
	@debug order moiety.atoms.species[order], moiety.atoms.coords[order]
	R_group = [length(order):-1:(length(order) - length(R_group_indices) + 1)...]

	# rebuild Atoms
	atoms = Atoms(moiety.atoms.species[order], moiety.atoms.coords[order]) # atoms are sorted by degree and un-tagged

	# build moiety with ordered atoms
	moiety = Crystal(name, box, atoms, charges)
    infer_bonds!(moiety, false)
	if length(R_group_indices) ≠ 0
		tag_R_group!(moiety, R_group_indices, df) # replace tags
	end

	return moiety # nodes are sorted by bond order, and R group is tagged w/ _
end

end # module
