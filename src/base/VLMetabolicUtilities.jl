# == PRIVATE METHODS BELOW HERE ======================================================================================= #
function _extract_species_from_phrase(phrase::String)::Array{String,1}

    # cut around the +'s
    species_array = Array{String,1}()
    species_substring_array = split(phrase, "+")

    # process each substring -
    for species_substring in species_substring_array
        
        # get the last value after we split for * -
        species_value = string(last(split(species_substring, "*")))
        
        # grab -
        push!(species_array, species_value)
    end

    # return -
    return species_array
end

function _extract_stochiometric_coefficient_from_phrase(phrase::String, species::String)::Float64

    # check: does this species occur in this phrase?
    if (occursin(species, phrase) == false)
        return 0.0
    end

    # ok: if we get here, then we have this species -
    # split around the +
    plus_split_product_array = split(phrase, "+") 

    # process each component of the phrase -
    for species_substring in plus_split_product_array
    
        # get the values after we split for * -
        tmp_value_array = string.(split(species_substring, "*"))

        # ok: the last element will always be the species. If we have two elements the first is the st coeff -
        if (last(tmp_value_array) == species)
            if (length(tmp_value_array) == 2)
                return parse(Float64, first(tmp_value_array))
            else
                return 1.0
            end
        end
    end

    # default returns 0.0 -
    return 0.0
end

function _process_left_phrase(phrase::String, species::String)::Float64

    # get coeff -
    stm_coeff = _extract_stochiometric_coefficient_from_phrase(phrase, species)

    # ok: if non-zero, then multiply by -1
    if (stm_coeff > 0)
        return -1.0 * stm_coeff
    end

    # default -
    return stm_coeff
end

function _process_right_phrase(phrase::String, species::String)::Float64
    return _extract_stochiometric_coefficient_from_phrase(phrase, species)
end

# == PRIVATE METHODS ABOVE HERE ======================================================================================= #

# == PUBLIC METHODS BELOW HERE ======================================================================================== #
function build_stoichiometric_matrix(reactionTable::DataFrame, speciesSymbolArray::Array{String,1})::Some

    try

        # initialize -
        reaction_id_array = build_reaction_id_array(reactionTable) |> check
        species_symbol_array = speciesSymbolArray
        number_of_reactions = length(reaction_id_array)
        number_of_species = length(species_symbol_array)
        stoichiometric_matrx = zeros(number_of_species, number_of_reactions)

        # main -
        for species_index = 1:number_of_species # rows
            
            # grab species symbol -
            species_symbol = species_symbol_array[species_index]
        
            # does this symbol appear in a reaction?
            for reaction_index = 1:number_of_reactions
                
                # get the L and R phrases -
                L = reactionTable[reaction_index,:forward]
                R = reactionTable[reaction_index,:reverse]

                # get the stcoeff from the L phrase -
                L_st_coeff = _process_left_phrase(L, species_symbol)
            
                # get the stcoeff from the R phrase -
                R_st_coeff = _process_right_phrase(R, species_symbol)
            
                # update the stoichiometric_matrx -
                stoichiometric_matrx[species_index,reaction_index] = (L_st_coeff + R_st_coeff)
            end
        end

        # return -
        return Some(stoichiometric_matrx)
    catch error
        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end

function build_stoichiometric_matrix(reactionTable::DataFrame)::Some

    try 

        # initialize -
        species_symbol_array = build_species_symbol_array(reactionTable) |> check
        stoichiometric_matrx = build_stoichiometric_matrix(reactionTable, species_symbol_array) |> check

        # return -
        return Some(stoichiometric_matrx)
    catch error
        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end

function build_flux_bounds_array(reactionTable::DataFrame; defaultFluxBoundValue::Float64=100.0)::Some

    try
        
        # initialize -
        reaction_id_array = build_reaction_id_array(reactionTable) |> check
        number_of_reactions = length(reaction_id_array)
        flux_bounds_array = Array{Float64,2}(undef, number_of_reactions, 2)

        # iterate through the reaction id's and determine if the reaction is reversible. If so, then build the bounds array 
        # using the default values -
        for (index, id_value) in enumerate(reaction_id_array)
            
            # is this reaction reversible?
            is_reversible = reactionTable[index,:reversibility]
            if (is_reversible == true)
                flux_bounds_array[index,1] = -1.0 * defaultFluxBoundValue
                flux_bounds_array[index,2] = defaultFluxBoundValue
            elseif (is_reversible == false)
                flux_bounds_array[index,1] = 0.0
                flux_bounds_array[index,2] = defaultFluxBoundValue
            end
        end

        # return -
        return Some(flux_bounds_array)
    catch error
        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end

function build_species_bounds_array(speciesTable::DataFrame)::Some

    try 

        # initialize -
        (number_of_species, _) = size(speciesTable)
        species_bounds_array = Array{Float64,2}(undef, number_of_species, 2)

        # main -
        for species_index = 1:number_of_species
            
            # get data for bounds array -
            lower_bound = speciesTable[species_index, :lower_bound]
            upper_bound = speciesTable[species_index, :upper_bound]

            # package -
            species_bounds_array[species_index,1] = lower_bound
            species_bounds_array[species_index,2] = upper_bound
        end

        # return -
        return Some(species_bounds_array)
    catch error
        
        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end

function build_species_table(reactionTable::DataFrame)::Some

    try

        # initialize -
        species_symbol_array = build_species_symbol_array(reactionTable) |> check
        df_species_table = DataFrame(symbol=String[], name=Union{String,Missing}[], lower_bound=Float64[], upper_bound=Float64[])

        # main -
        for species_symbol in species_symbol_array
            
            # create a tuple w/the data row -
            data_row = (species_symbol, missing, 0.0, 0.0)
            
            # add to the df -
            push!(df_species_table, data_row)
        end

        # return -
        return Some(df_species_table)
    catch error

        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end

function build_species_symbol_array(reactionTable::DataFrame)::Some

    try
        
        # initialize -
        species_symbol_array = Array{String,1}()
        reaction_phrase_array = Array{String,1}()
        (number_of_reactions, _) = size(reactionTable)
        
        # populate the array of reaction phrases -
        for reaction_index = 1:number_of_reactions
            
            # get the L and R phrases -
            L = reactionTable[reaction_index,:forward]
            R = reactionTable[reaction_index,:reverse]

            # capture -
            push!(reaction_phrase_array, L)
            push!(reaction_phrase_array, R)
        end

        # ok, so now lets populate the tmp species set -
        for reaction_phrase in reaction_phrase_array
            
            # get species sub array -
            species_in_reaction_phrase = _extract_species_from_phrase(reaction_phrase)

            # push into set -
            for tmp_species_symbol in species_in_reaction_phrase
                if (in(tmp_species_symbol, species_symbol_array) == false)
                    push!(species_symbol_array, tmp_species_symbol)
                end
            end
        end

        # return -
        return Some(species_symbol_array)
    catch error
        
        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end

function build_objective_coefficient_array(reactionTable::DataFrame)::Some

    try

        # initialize -
        (number_of_reactions, _) = size(reactionTable)

        # default behavior -
        obj_coeff_array = zeros(number_of_reactions)
        
        # return -
        return Some(obj_coeff_array)

    catch error

        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end

function build_reaction_id_array(reactionTable::DataFrame)::Some

    try
        
        # get the id col -
        id_col = reactionTable[!,:id]
    
        # return -
        return Some(id_col)
    catch error
        
        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end

function build_transport_reaction_table(species_pair_array::Array{Pair{String,String},1})::Some

    try

        # initailize -
        id_array = String[]
        forward_reaction_string = String[]
        reverse_reaction_string = String[]
        reversibility_array = Bool[]
        ec_number_array = Union{Missing,String}[]
        lower_bound_array = Union{Missing,Float64}[]
        upper_bound_array = Union{Missing,Float64}[]

        # ok, so we have an array of pairs - assume all transport reactions are *reversible* -
        for transport_pair in species_pair_array
            
            # get the a => b
            first_species = transport_pair.first
            second_species = transport_pair.second

            # write the record -
            push!(id_array, "$(first_species)_exchange")
            push!(forward_reaction_string, first_species)
            push!(reverse_reaction_string, second_species)
            push!(reversibility_array, true)
            push!(ec_number_array, missing)
            push!(lower_bound_array, -Inf)
            push!(upper_bound_array, Inf)
        end

        # package into DataFrame -
        reaction_dataframe = DataFrame(id=id_array, forward=forward_reaction_string, reverse=reverse_reaction_string, reversibility=reversibility_array, 
            ec=ec_number_array, lower_bound=lower_bound_array, upper_bound=upper_bound_array)

        # return -
        return Some(reaction_dataframe)
    catch error

        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end

function build_atom_matrix(formula_array::Array{String,1}; 
        elements=["C","H","N","O","P","S"])::Some
    
    try

        # initialize -
        number_of_elements = length(elements)
        number_of_compounds = length(formula_array)
        atom_matrix = Array{Number,2}(undef,number_of_elements,number_of_compounds)

        # ok, so let's process this list of formulas -
        for (formula_index,formula) in enumerate(formula_array)

            # parse -
            atom_dictionary = parse_molecular_formula_string(formula) |> check

            # fill in the elements -
            for (element_index,element) in enumerate(elements)
                
                # check: do we have this element?
                get!(atom_dictionary, element, 0)
                
                # go -
                atom_matrix[element_index, formula_index] = atom_dictionary[element]
            end
        end

        # return -
        return Some(atom_matrix)
    catch error

        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end

function parse_molecular_formula_string(formula::String)::Some
    
    try

        # test -
        atom_dictionary = Dict();
        local_array = Array{Char,1}()

        # turn string into char array -
        formula_char_array = collect(formula);

        # add extra 1 if last char is a letter -
        if (isnumeric(last(formula_char_array))== false)
            push!(formula_char_array,'1')
        end

        # read from the bottom -
        reverse!(formula_char_array)

        # how many chars do we have?
        while (isempty(formula_char_array) == false)
            
            # clean out the array from the last pass -
            empty!(local_array)

            # grab the next value -
            next_value = pop!(formula_char_array)
            if (isnumeric(next_value) == false)
                
                # we have an element -> read until I hit another element -
                is_ok_to_loop = true
                while (is_ok_to_loop)

                    if (isempty(formula_char_array) == true)
                        break
                    end

                    read_one_ahead = pop!(formula_char_array)
                    if (isnumeric(read_one_ahead) == true)
                        push!(local_array,read_one_ahead)
                    else

                        # ok: so if we get here - then we read the next char, but it was a 
                        # letter (element) - so we need to puch it back on the stack
                        push!(formula_char_array,read_one_ahead)

                        # reverse -
                        is_ok_to_loop = false
                    end
                end

                # we need to turn local array into a string -
                buffer = ""
                [buffer*=string(x) for x in local_array]
                atom_dictionary[string(next_value)] = parse(Int64,buffer)
            end
        end
        
        # return -
        return Some(atom_dictionary)

    catch error

        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end

# == PUBLIC METHODS ABOVE HERE ======================================================================================== #