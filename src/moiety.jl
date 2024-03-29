"""
Returns bonding rules involving R-group-tagged atoms
"""
function tagged_bonding_rules()::Array{BondingRule}
    newrules = []
    for rule in rc[:bonding_rules]
        if rule.species_i != :*
            push!(
                newrules,
                BondingRule(Symbol("$(rule.species_i)!"), rule.species_j, rule.max_dist)
            )
            push!(
                newrules,
                BondingRule(Symbol("$(rule.species_j)!"), rule.species_i, rule.max_dist)
            )
            push!(
                newrules,
                BondingRule(
                    Symbol("$(rule.species_i)!"),
                    Symbol("$(rule.species_j)!"),
                    rule.max_dist
                )
            )
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
    for i in r
        xtal.atoms.species[i] = Symbol(split("$(xtal.atoms.species[i])", rc[:r_tag])[1])
    end
end

"""
Returns a copy of a crystal w/ R group atoms deleted
"""
function subtract_r_group(xtal::Crystal)::Crystal
    not_r = [i for i in eachindex(xtal.atoms.species) if !(i ∈ r_group_indices(xtal))]
    coords = xtal.atoms.coords[not_r]
    species = xtal.atoms.species[not_r]
    return Crystal("no_r_$(xtal.name)", xtal.box, Atoms(species, coords), xtal.charges)
end

## moiety import function (exposed)
@doc raw"""

    q = moiety(xyz_filename)


Generates a moiety (`Crystal`) from an .xyz file found in `rc[:paths][:moieties]`.

Use `set_path_to_data` or set `rc[:paths][:moieties]` to change the path from which the XYZ file is read.

Atoms appended with '!' are tagged for replacement via `substructure_replace`.

Bonds are inferred automatically via `infer_bonds!`.

# Arguments
- `xyz_filename::Union{String,Nothing}` the moiety input file name, an `.xyz` file; if set to `nothing` the moiety is the null set.
- `bonding_rules::Union{Vector{BondingRule},Nothing}` (optional) a list of rules to use for inferring the bonding network of the atoms loaded from the XYZ file. If set to `nothing`, the default rules are used.
- `presort::Bool` whether to sort the atoms by bonding order for structure search efficiency. Set `false` to skip pre-sorting and maintain indexing order with source file. Does not apply to !-tagged atoms, which will still be moved to the end of the list.
"""
function moiety(
    name::Union{String, Nothing};
    bonding_rules::Union{Vector{BondingRule}, Nothing}=nothing,
    presort::Bool=true
)::Crystal
    # make box (arbitrary unit cube)
    box = unit_cube()
    # handle deletion option (replace-with-nothing)
    if !isnothing(name)
        xf = Frac(read_xyz("$(rc[:paths][:moieties])/$name"), box)
    else
        name = "nothing"
        xf = Atoms{Frac}(0)
    end
    # generate Crystal from moiety XYZ coords
    charges = Charges{Frac}(0)
    moiety = Crystal(name, box, xf, charges)
    # ID R group
    R_group_indices = r_group_indices(moiety)
    # handle custom vs. default bonding rules
    if isnothing(bonding_rules)
        infer_bonds!(moiety, false)
    else
        infer_bonds!(moiety, false; bonding_rules=bonding_rules)
    end
    # sort by node degree
    sp = sortperm(degree(moiety.bonds); rev=true)
    order = presort ? sp : eachindex(sp)
    # ordered atoms
    if length(R_group_indices) > 0
        order_wo_R = order[[i for i in eachindex(order) if !(order[i] ∈ R_group_indices)]]
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
