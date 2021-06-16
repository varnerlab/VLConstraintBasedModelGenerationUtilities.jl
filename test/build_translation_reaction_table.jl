using VLConstraintBasedModelGenerationUtilities
using DataFrames

# setup path to sequence file -
path_to_protein_sequence_file = "/Users/jeffreyvarner/Desktop/julia_work/VLConstraintBasedModelGenerationUtilities.jl/test/data/P-MZ373340.fasta"

# load -
protein_table = build_protein_table(path_to_protein_sequence_file) |> check
seq = protein_table[1,:protein_sequence]
seq_name = protein_table[1,:id]

# build the translation table -
translation_dictionary = Dict{String,DataFrame}()
(number_of_proteins, _) = size(protein_table)
for protein_index = 1:number_of_proteins
    
    local_seq = protein_table[protein_index,:protein_sequence]
    local_seq_name = protein_table[protein_index,:id]
    translation_table = build_translation_reaction_table(local_seq_name, local_seq) |> check
    
    # grab -
    translation_dictionary[local_seq_name] = translation_table
end 
