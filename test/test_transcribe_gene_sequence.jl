using VLConstraintBasedModelGenerationUtilities
using BioSymbols

# setup path to sequence file -
path_to_gene_sequence_file = "/Users/jeffreyvarner/Desktop/julia_work/VLConstraintBasedModelGenerationUtilities.jl/test/data/G-MZ373340.fasta"

# load -
gene_table = build_gene_table(path_to_gene_sequence_file) |> check

# get the sequence -
gene_seq = gene_table[!,:gene_sequence][1]

# setup local complementOp function -
function complementOp(nucleotide::BioSymbols.DNA)::BioSymbols.RNA

    if nucleotide == DNA_T
        return RNA_U
    elseif nucleotide == DNA_A
        return RNA_A
    elseif nucleotide == DNA_C
        return RNA_C
    elseif nucleotide == DNA_G
        return RNA_G
    else
        return RNA_N
    end
end

# table -
transcription_table = transcribe_sequence(gene_seq; complementOperation=complementOp) |> check