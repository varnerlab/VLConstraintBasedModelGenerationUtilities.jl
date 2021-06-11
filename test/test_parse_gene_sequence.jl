using VLConstraintBasedModelGenerationUtilities

# setup path to sequence file -
path_to_gene_sequence_file = "/Users/jeffreyvarner/Desktop/julia_work/VLConstraintBasedModelGenerationUtilities.jl/test/data/test.fasta"

# load -
result = build_gene_table(path_to_gene_sequence_file) |> check