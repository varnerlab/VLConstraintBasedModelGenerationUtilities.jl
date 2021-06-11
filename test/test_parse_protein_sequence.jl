using VLConstraintBasedModelGenerationUtilities

# setup path to sequence file -
path_to_protein_sequence_file = "/Users/jeffreyvarner/Desktop/julia_work/VLConstraintBasedModelGenerationUtilities.jl/test/data/P-MZ373340.fasta"

# load -
result = build_protein_table(path_to_protein_sequence_file) |> check