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
export build_transport_reaction_table
export transcribe
export translate
export count

# metabolic methods -
export build_reaction_id_array
export build_flux_bounds_array
export build_species_symbol_array

# export complement function -
export !

end # module
