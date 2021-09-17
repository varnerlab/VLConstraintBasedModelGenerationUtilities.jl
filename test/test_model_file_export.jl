using VLConstraintBasedModelGenerationUtilities

# setup path to protein sequence file -
path_to_vff_file = "/Users/jeffreyvarner/Desktop/julia_work/VLConstraintBasedModelGenerationUtilities.jl/test/data/Test.vff"
path_to_system_model_file = "/Users/jeffreyvarner/Desktop/julia_work/VLConstraintBasedModelGenerationUtilities.jl/test/data/Test.bson"

# let's build the reaction table -
metabolic_reaction_table = build_metabolic_reaction_table(path_to_vff_file) |> check

# build the stm -
stm = build_stoichiometric_matrix(metabolic_reaction_table) |> check

# build the species bounds array -
species_table = build_species_table(metabolic_reaction_table) |> check
species_bounds_array = build_species_bounds_array(species_table) |> check

# build the flux bounds array -
flux_bounds_array = build_flux_bounds_array(metabolic_reaction_table) |> check

# write the system model file -
result = write_system_model_file(path=path_to_system_model_file, stoichiometric_matrix=stm, 
    flux_bounds_array=flux_bounds_array, species_bounds_array=species_bounds_array)

