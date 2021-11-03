module VLConstraintBasedModelGenerationUtilities

# Include my files -
include("Include.jl")

# export methods and types -
export check
export build_gene_table
export build_protein_table
export build_metabolic_reaction_table
export build_transcription_reaction_table
export build_translation_reaction_table
export build_transport_reaction_table
export build_txtl_transport_reaction_table
export transcribe
export translate
export count

# metabolic methods -
export build_reaction_id_array
export build_flux_bounds_array
export build_species_symbol_array
export build_species_bounds_array
export build_species_table
export build_species_symbol_table
export build_stoichiometric_matrix
export build_atom_matrix
export build_objective_coefficient_array
export parse_molecular_formula_string

# file io methods -
export write_system_model_file

# private: but has a unit test
# export extract_stochiometric_coefficient_from_phrase

# export complement function -
export !

end # module
