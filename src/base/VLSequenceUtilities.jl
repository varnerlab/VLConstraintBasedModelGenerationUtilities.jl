function build_transcription_reaction_table(table::DataFrame; polymeraseSymbol::Symbol=:RNAP, 
    logger::Union{Nothing,SimpleLogger}=nothing)::VLResult

    try
    catch error
        return VLResult(error)
    end
end

function build_translation_reaction_table(table::DataFrame; ribosomeSymbol::Symbol=:RIBOSOME, 
    logger::Union{Nothing,SimpleLogger}=nothing)::VLResult

    try
    catch error
        return VLResult(error)
    end
end

function transcribe_sequence(sequence::BioSequences.LongSequence; 
    logger::Union{Nothing,SimpleLogger}=nothing)

    try
    catch error
        return VLResult(error)
    end
end

function translate_sequence(sequence::BioSequences.LongSequence; 
    logger::Union{Nothing,SimpleLogger}=nothing)


    try
    catch error
        return VLResult(error)
    end
end

function transcribe_sequence(table::DataFrame; 
    logger::Union{Nothing,SimpleLogger}=nothing)

    try

        # initailize -
        transcription_dictionary = Dict{String,Any}()

        # get the length of the table -
        (number_of_rows, _) = size(table)
        for row_index = 1:number_of_rows
            
            # get id -
            sequence_id = table[row_index, :id]
            sequence = table[row_index, :sequence]

            # transcribe -
            result = transcribe_sequence(sequence; logger=logger) |> check

            # package -
            transcription_dictionary[sequence_id] = result
        end

        # return -
        return VLResult(transcription_dictionary)
    catch error
        return VLResult(error)
    end
end

function translate_sequence(table::DataFrame; 
    logger::Union{Nothing,SimpleLogger}=nothing)


    try
    catch error
        return VLResult(error)
    end
end