using VLConstraintBasedModelGenerationUtilities
using DataFrames

# setup path to sequence file -
path_to_gene_sequence_file = "/Users/jeffreyvarner/Desktop/julia_work/VLConstraintBasedModelGenerationUtilities.jl/test/data/G-MZ373340.fasta"

# load -
gene_table = build_gene_table(path_to_gene_sequence_file) |> check

# get the sequence -
gene_seq = gene_table[!,:gene_sequence][1]
gene_seq_name = gene_table[!,:id][1]

# transcription table -
transcription_table = build_transcription_reaction_table(gene_seq_name, gene_seq) |> check

# setup path to protein sequence file -
path_to_protein_sequence_file = "/Users/jeffreyvarner/Desktop/julia_work/VLConstraintBasedModelGenerationUtilities.jl/test/data/P-MZ373340.fasta"

# load -
protein_table = build_protein_table(path_to_protein_sequence_file) |> check

# build the translation table -
translation_dictionary = Dict{String,DataFrame}()
translation_array = Array{DataFrame,1}()
(number_of_proteins, _) = size(protein_table)
for protein_index = 1:number_of_proteins
    
    local_seq = protein_table[protein_index,:protein_sequence]
    local_seq_name = protein_table[protein_index,:id]
    translation_table = build_translation_reaction_table(local_seq_name, local_seq) |> check
    
    # grab -
    # translation_dictionary[local_seq_name] = translation_table
    push!(translation_array, translation_table)
end

# add together -
master_translation_table = reduce(vcat, translation_array)

# lastly - let's build the transport reactions -
transport_reaction_table = build_transport_reaction_table() |> check

# append reactions into master reaction table -
master_reaction_table = reduce(vcat, [transcription_table, master_translation_table, transport_reaction_table])

