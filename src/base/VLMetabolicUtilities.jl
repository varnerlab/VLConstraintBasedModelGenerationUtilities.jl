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
function build_stoichiometric_matrix(reactionTable::DataFrame)::VLResult

    try 

        # initialize -
        reaction_id_array = build_reaction_id_array(reactionTable) |> check
        species_symbol_array = build_species_symbol_array(reactionTable) |> check
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
        return VLResult(stoichiometric_matrx)
    catch error
        return VLResult(error)
    end
end

function build_flux_bounds_array(reactionTable::DataFrame; defaultFluxBoundValue::Float64=100.0)::VLResult

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
        return VLResult(flux_bounds_array)
    catch error
        return VLResult(error)
    end
end

function build_species_bounds_array()::VLResult

    try 
    catch error
        return VLResult(error)
    end
end

function build_species_symbol_array(reactionTable::DataFrame)::VLResult

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
        return VLResult(species_symbol_array)
    catch error
        return VLResult(error)
    end
end

function build_reaction_id_array(reactionTable::DataFrame)::VLResult

    try
        
        # get the id col -
        id_col = reactionTable[!,:id]
    
        # return -
        return VLResult(id_col)
    catch error
        return VLResult(error)
    end
end
# == PUBLIC METHODS ABOVE HERE ======================================================================================== #