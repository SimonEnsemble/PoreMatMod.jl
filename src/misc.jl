# Retuns the geometric center of an Array, Frac/Atoms object, or Crystal.
function geometric_center(xf::Array{Float64,2})::Array{Float64}
    return sum(xf, dims=2)[:] / size(xf, 2)
end

geometric_center(coords::Frac)::Array{Float64} = geometric_center(coords.xf)

geometric_center(atoms::Atoms)::Array{Float64} = geometric_center(atoms.coords)

geometric_center(xtal::Crystal)::Array{Float64} = geometric_center(xtal.atoms)

# Helper for making .xyz's
write_xyz(xtal::Crystal, name::String) = Xtals.write_xyz(Cart(xtal.atoms, xtal.box), name)


# Translates all atoms in xtal such that xtal[1] is in its original position
# and the rest of xtal is in its nearest-image position relative to xtal[1]
function adjust_for_pb!(xf::Matrix{Float64}; xtal_name::String="")
    # record position vector of xtal[1]
    origin_offset = deepcopy(xf[:, 1])
    # loop over atom indices and correct coordinates
    for i in 1:size(xf, 2)
        # move atoms near the origin for nearest-image calculation
        dxf = xf[:, i] .- origin_offset
        # nearest_image! expects points to be within same or adjacent unit cells
        @assert all(abs.(dxf) .< 2) "Invalid xf coords $xtal_name"
        # resolve periodic boundaries (other vectors unchanged)
        nearest_image!(dxf)
        # return atoms to their [nearest-image] original positions
        xf[:, i] = dxf .+ origin_offset
    end
end

adjust_for_pb!(xtal::Crystal) = adjust_for_pb!(xtal.atoms.coords.xf, xtal_name=xtal.name)


# returns an Array containing the indices
function idx_filter(xtal::Crystal, subset::Array{Int})::Array{Int,1}
    return [i for i in 1:xtal.atoms.n if !(i ∈ subset)]
end
