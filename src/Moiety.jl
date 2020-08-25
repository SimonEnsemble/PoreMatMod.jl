PATH_TO_MOIETIES = joinpath(pwd(), "data/moieties")
R_GROUP_TAG = '!'


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
	return R
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
	@debug 1 moiety.atoms.species
	# ID R group
	R_group_indices = filter_R_group(moiety)
	@debug 2 R_group_indices
	# sort by node degree
	bondingrules = new_bonding_rules()
	infer_bonds!(moiety, false, bondingrules)
	order = sortperm(degree(moiety.bonds), rev=true)
	@debug 3 order
	# ordered atoms
	order_wo_R = length(R_group_indices) > 0 ? order[1:end .!= R_group_indices] : order # BUG this is wrong. it is removing the element at an index instead of a value from an array
	@debug 4 order_wo_R
	# append R-group to the end
	order = vcat(order_wo_R, R_group_indices)
	@debug 5 order
	# rebuild Atoms
	atoms = Atoms(moiety.atoms.species[order], moiety.atoms.coords[order])
	@debug 6 atoms
	# nodes are sorted by bond order, and R group is moved to end & tagged w/ !
	moiety = Crystal(name, box, atoms, charges)
	@debug 7 moiety
	infer_bonds!(moiety, false, bondingrules)
	return moiety
end
