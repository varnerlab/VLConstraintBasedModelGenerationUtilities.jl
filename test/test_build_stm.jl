using VLConstraintBasedModelGenerationUtilities

# setup path to protein sequence file -
path_to_vff_file = "/Users/jeffreyvarner/Desktop/julia_work/VLConstraintBasedModelGenerationUtilities.jl/test/data/Test.vff"

# let's build the reaction table -
metabolic_reaction_table = build_metabolic_reaction_table(path_to_vff_file) |> check

# get list of species -
list_of_species = build_species_symbol_array(metabolic_reaction_table) |> check

# build the stm -
stm = build_stoichiometric_matrix(metabolic_reaction_table) |> check