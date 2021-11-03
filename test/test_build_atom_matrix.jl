using VLConstraintBasedModelGenerationUtilities

# build an array of compounds -
formula_array = [
    "C6H12O6"        ;      # 1 glucose
    "C6H13O9P"       ;      # 2 g6p
    "C10H16N5O13P3"  ;      # 3 atp
    "C10H15N5O10P2"  ;      # 4 adp
];

# test -
compound_dict = parse_molecular_formula_string(formula_array[4]) |> check

# compute the atom matrix -
atom_matrix = build_atom_matrix(formula_array) |> check