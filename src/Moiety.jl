PATH_TO_MOIETIES = joinpath(pwd(), "data/moieties")
R_GROUP_TAG = '!'

"""
Returns bonding rules including R-group-tagged atom copies
"""
function new_bonding_rules()::Array{BondingRule}
    bondingrules = Xtals.bondingrules()
    push!(bondingrules, BondingRule(:C, :*, 0.4, 1.9))
    newrules = []
    for i in 1:length(bondingrules)
        if bondingrules[i].species_i != :*
            push!(newrules, BondingRule(
                Symbol("$(bondingrules[i].species_i)!"),
                bondingrules[i].species_j,
                bondingrules[i].min_dist,
                bondingrules[i].max_dist))
            push!(newrules, BondingRule(
                Symbol("$(bondingrules[i].species_j)!"),
                bondingrules[i].species_i,
                bondingrules[i].min_dist,
                bondingrules[i].max_dist))
            push!(newrules, BondingRule(
                Symbol("$(bondingrules[i].species_i)!"),
                Symbol("$(bondingrules[i].species_j)!"),
                bondingrules[i].min_dist,
                bondingrules[i].max_dist))
        end
    end
    return unique(vcat(bondingrules, newrules))
end


"""
Returns R group indices (whichever atoms have species symbols appended by '!')
"""
function r_group_indices(xtal::Crystal)::Array{Int}
    @debug "Filtering R group" R_GROUP_TAG
    R = []
    for (idx, label) in enumerate(xtal.atoms.species) # loop over crystal atoms to find tags
        # if String representation of label Symbol ends in !, atom is in R
        tokens = split("$label", R_GROUP_TAG)
        if length(tokens) == 2 && tokens[2] == "" # other ! in symbol not tolerated.
            push!(R, idx)
        end
    end
    return R
end


"""
Un-tags R group atoms (removes '!' suffix)
"""
function untag_r_group!(xtal::Crystal)
    @debug "Un-tagging R group in $(xtal.name)" R_GROUP_TAG
    r = r_group_indices(xtal) # get indices of R group
    for i ∈ r
        xtal.atoms.species[i] = Symbol(split("$(xtal.atoms.species[i])", R_GROUP_TAG)[1])
    end
end


"""
Returns a copy of a crystal w/ R group atoms deleted
"""
function subtract_r_group(xtal::Crystal)::Crystal
    not_r = [i for i in 1:length(xtal.atoms.species) if !(i ∈ r_group_indices(xtal))]
    coords = xtal.atoms.coords[not_r]
    species = xtal.atoms.species[not_r]
    return Crystal("no_r_$(xtal.name)", xtal.box, Atoms(species, coords), xtal.charges)
end


## moiety import function (exposed)
@doc raw"""
    moiety(name)

Generates a moiety (`Crystal`) from an .xyz file found in `PATH_TO_MOIETIES`

Use `set_path_to_moieties` or `set_path_to_data` to change the input path.

Atoms appended with '!' are tagged for replacement via `find_replace`.

Bonds are inferred within the local unit cell only (no bonds across periodic boundaries).

# Arguments
- `name::String` the moiety name (input file name without the .xyz extension).
"""
function moiety(name::Union{String,Nothing})::Crystal
    @debug "Getting moiety: $name"
    # generate Crystal from moiety XYZ coords
    box = unit_cube()
    if !isnothing(name)
        fx = Frac(read_xyz(joinpath(pwd(), "$PATH_TO_MOIETIES/$name.xyz")), box)
    else
        name = "nothing"
        fx = Atoms{Frac}(0)
    end
    charges = Charges{Frac}(0)
    moiety = Crystal(name, box, fx, charges)
    # ID R group
    R_group_indices = r_group_indices(moiety)
    # sort by node degree
    bondingrules = new_bonding_rules()
    infer_bonds!(moiety, false, bonding_rules=bondingrules)
    order = sortperm(degree(moiety.bonds), rev=true)
    # ordered atoms
    if length(R_group_indices) > 0
        order_wo_R = order[
            [i for i in 1:length(order) if !(order[i] ∈ R_group_indices)]]
    else
        order_wo_R = order
    end
    # append R-group to the end
    order = vcat(order_wo_R, R_group_indices)
    # rebuild Atoms
    atoms = Atoms(moiety.atoms.species[order], moiety.atoms.coords[order])
    # nodes are sorted by bond order, and R group is moved to end & tagged w/ !
    moiety = Crystal(name, box, atoms, charges)
    infer_bonds!(moiety, false, bonding_rules=bondingrules)
    return moiety
end
