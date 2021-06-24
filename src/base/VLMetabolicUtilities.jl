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

# == PRIVATE METHODS ABOVE HERE ======================================================================================= #

# == PUBLIC METHODS BELOW HERE ======================================================================================== #
function build_stoichiometric_matrix(reactionTable::DataFrame)::VLResult

    try 
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