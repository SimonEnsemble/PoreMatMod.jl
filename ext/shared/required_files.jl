required_files = Dict(
        :crystals => [
            "IRMOF-1.cif",
            "SIFSIX-2-Cu-i.cif",
            "IRMOF-1_noH.cif",
            "UiO-66.cif",
            "NiPyC_fragment_trouble.cif"
        ],
        :moieties => [
            "2-!-p-phenylene.xyz",
            "2-acetylamido-p-phenylene.xyz",
            "1,4-C-phenylene_noH.xyz",
            "1,4-C-phenylene.xyz",
            "4-pyridyl.xyz",
            "acetylene.xyz",
            "BDC.xyz",
            "disordered_ligand!.xyz",
            "formate_caps.xyz",
            "SBU.xyz"
        ]
)

function check_example_data()
    # make sure directories are present and the right files for the examples
    for file_type in [:moieties, :crystals]
        # make sure directories exist
        if !isdir(rc[:paths][file_type])
            @warn "$(rc[:paths][file_type]) directory not present; creating it."
            mkpath(rc[:paths][file_type])
        end
        for required_file in required_files[file_type]
            where_it_shld_be = joinpath(rc[:paths][file_type], required_file)
            if !isfile(where_it_shld_be)
                @warn "$where_it_shld_be not present; copying it from src."
                where_it_is = normpath(
                    joinpath(
                        @__DIR__,
                        "..",
                        "..",
                        "examples",
                        "data",
                        String(file_type),
                        required_file
                    )
                )
                cp(where_it_is, where_it_shld_be)
            end
        end
    end
end
