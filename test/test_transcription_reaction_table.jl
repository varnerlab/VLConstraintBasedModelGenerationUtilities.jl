using VLConstraintBasedModelGenerationUtilities

# setup path to sequence file -
path_to_gene_sequence_file = "/Users/jeffreyvarner/Desktop/julia_work/VLConstraintBasedModelGenerationUtilities.jl/test/data/G-MZ373340.fasta"

# load -
gene_table = build_gene_table(path_to_gene_sequence_file) |> check

# get the sequence -
gene_seq = gene_table[!,:gene_sequence][1]

# table -
transcription_table = build_transcription_reaction_table("test_gene", gene_seq) |> check