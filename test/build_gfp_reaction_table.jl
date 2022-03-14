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