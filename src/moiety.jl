"""
Returns bonding rules involving R-group-tagged atoms
"""
function tagged_bonding_rules()::Array{BondingRule}
    newrules = []
    for rule ∈ rc[:bonding_rules]
        if rule.species_i != :*
            push!(newrules, BondingRule(Symbol("$(rule.species_i)!"), rule.species_j, rule.max_dist))
            push!(newrules, BondingRule(Symbol("$(rule.species_j)!"), rule.species_i, rule.max_dist))
            push!(newrules, BondingRule(Symbol("$(rule.species_i)!"), Symbol("$(rule.species_j)!"), rule.max_dist))
        end
    end
    return newrules
end


"""
Returns R group indices (whichever atoms have species symbols appended by '!')
"""
function r_group_indices(xtal::Crystal)::Array{Int}
    R = []
    for (idx, label) in enumerate(xtal.atoms.species) # loop over crystal atoms to find tags
        # if String representation of label Symbol ends in !, atom is in R
        tokens = split("$label", rc[:r_tag])
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
    r = r_group_indices(xtal) # get indices of R group
    for i ∈ r
        xtal.atoms.species[i] = Symbol(split("$(xtal.atoms.species[i])", rc[:r_tag])[1])
    end
end


"""
Returns a copy of a crystal w/ R group atoms deleted
"""
function subtract_r_group(xtal::Crystal)::Crystal
    not_r = [i for i ∈ 1:length(xtal.atoms.species) if !(i ∈ r_group_indices(xtal))]
    coords = xtal.atoms.coords[not_r]
    species = xtal.atoms.species[not_r]
    return Crystal("no_r_$(xtal.name)", xtal.box, Atoms(species, coords), xtal.charges)
end


## moiety import function (exposed)
@doc raw"""
```julia
moiety(xyz_filename)
```

Generates a moiety (`Crystal`) from an .xyz file found in `rc[:paths][:moieties]`.

Use `set_path_to_data` or set `rc[:paths][:moieties]` to change the path from which the XYZ file is read.

Atoms appended with '!' are tagged for replacement via `substructure_replace`.

Bonds are inferred automatically via `infer_bonds!`.

# Arguments
- `xyz_filename::String` the moiety input file name, an `.xyz` file.
"""
function moiety(name::Union{String,Nothing})::Crystal
    # generate Crystal from moiety XYZ coords
    box = unit_cube()
    if !isnothing(name)
        fx = Frac(read_xyz("$(rc[:paths][:moieties])/$name"), box)
    else
        name = "nothing"
        fx = Atoms{Frac}(0)
    end
    charges = Charges{Frac}(0)
    moiety = Crystal(name, box, fx, charges)
    # ID R group
    R_group_indices = r_group_indices(moiety)
    # sort by node degree
    infer_bonds!(moiety, false)
    order = sortperm(degree(moiety.bonds), rev=true)
    # ordered atoms
    if length(R_group_indices) > 0
        order_wo_R = order[
            [i for i ∈ 1:length(order) if !(order[i] ∈ R_group_indices)]]
    else
        order_wo_R = order
    end
    # append R-group to the end
    order = vcat(order_wo_R, R_group_indices)
    # rebuild Atoms
    atoms = Atoms(moiety.atoms.species[order], moiety.atoms.coords[order])
    # nodes are sorted by bond order, and R group is moved to end & tagged w/ !
    moiety = Crystal(name, box, atoms, charges)
    infer_bonds!(moiety, false)
    return moiety
end
