module VLConstraintBasedModelGenerationUtilities

# Include my files -
include("Include.jl")

# export methods and types -
export VLResult
export check
export build_gene_table
export build_protein_table
export build_metabolic_reaction_table
export build_transcription_reaction_table
export build_translation_reaction_table
export transcribe_sequence
export translate_sequence

# export complement function -
export !

end # module
