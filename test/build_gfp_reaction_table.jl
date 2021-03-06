using VLConstraintBasedModelGenerationUtilities
using DataFrames

# setup path to sequence file -
path_to_gene_sequence_file = "/Users/jeffreyvarner/Desktop/julia_work/VLConstraintBasedModelGenerationUtilities.jl/test/data/G-GFP.fasta"

# load -
gene_table = build_gene_table(path_to_gene_sequence_file) |> check

# get the sequence -
gene_seq = gene_table[!, :gene_sequence][1]
gene_seq_name = gene_table[!, :id][1]

# transcription table -
transcription_table = build_transcription_reaction_table(gene_seq_name, gene_seq; polymeraseSymbol = :RNAP) |> check

# setup path to protein sequence file -
path_to_protein_sequence_file = "/Users/jeffreyvarner/Desktop/julia_work/VLConstraintBasedModelGenerationUtilities.jl/test/data/P-GFP.fasta"

# load -
protein_table = build_protein_table(path_to_protein_sequence_file) |> check

local_seq = protein_table[1, :protein_sequence]
local_seq_name = protein_table[1, :id]
translation_table = build_translation_reaction_table(local_seq_name, local_seq) |> check

# append reactions into master reaction table -
master_reaction_table = reduce(vcat, [transcription_table, translation_table])

# get list of species -
list_of_species = build_species_symbol_array(master_reaction_table; left = :forward_reaction, right = :reverse_reaction) |> check

# build the stoichiometric_matrix -
stm = build_stoichiometric_matrix(master_reaction_table; left = :forward_reaction, right = :reverse_reaction) |> check

# build the species bounds array -
species_table = build_species_table(master_reaction_table; left = :forward_reaction, right = :reverse_reaction) |> check
species_bounds_array = build_species_bounds_array(species_table) |> check

# build the flux bounds array -
flux_bounds_array = build_flux_bounds_array(master_reaction_table) |> check

# write the system model file -
path_to_system_model_file = "/Users/jeffreyvarner/Desktop/julia_work/VLConstraintBasedModelGenerationUtilities.jl/test/data/GFP.bson"
result = write_system_model_file(path = path_to_system_model_file, stoichiometric_matrix = stm,
    flux_bounds_array = flux_bounds_array, species_bounds_array = species_bounds_array, reaction_table = master_reaction_table, species_table = species_table)