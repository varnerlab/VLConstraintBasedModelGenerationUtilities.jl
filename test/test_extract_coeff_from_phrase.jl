using VLConstraintBasedModelGenerationUtilities

# test phrase -
reaction_phrase = "3*M_a_c+8*M_b_c+M_d_c"

# test:
species = "M_gtp_c"
stm_coeff = extract_stochiometric_coefficient_from_phrase(reaction_phrase, species)