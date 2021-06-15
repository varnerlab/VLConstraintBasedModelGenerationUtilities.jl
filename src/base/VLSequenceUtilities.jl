function build_transcription_reaction_table(table::DataFrame; polymeraseSymbol::Symbol=:RNAP, 
    ecnumber::String="EC 2.7.7.6", logger::Union{Nothing,SimpleLogger}=nothing)::VLResult

    try

        # initialize -
        reaction_table = DataFrame(id=String[], forward=String[], reverse=String[], )

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

function transcribe_sequence(sequence::BioSequences.LongSequence, complementOperation::Function=!, 
    logger::Union{Nothing,SimpleLogger}=nothing)::VLResult

    try

        # initialize -
        new_sequence_array = Array{BioSymbol,1}()

        # ok: let's iterate the sequence, call the complementOperation function for each nt -
        for nt in sequence

            # call the complementOperation function to get the complementary nucleotide -
            complementary_nt = complementOperation(nt)

            # push -
            push!(new_sequence_array, complementary_nt)
        end

        # create new sequence -
        new_seq_from_array = LongSequence{DNAAlphabet{4}}(new_sequence_array)

        # return -
        return VLResult(new_seq_from_array)
    catch error
        return VLResult(error)
    end
end

function transcribe_sequence(table::DataFrame, complementOperation::Function=!; 
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
            result = transcribe_sequence(sequence, complementOperation; logger=logger) |> check

            # package -
            transcription_dictionary[sequence_id] = result
        end

        # return -
        return VLResult(transcription_dictionary)
    catch error
        return VLResult(error)
    end
end

function translate_sequence(sequence::BioSequences.LongSequence, complementOperation::Function;
    logger::Union{Nothing,SimpleLogger}=nothing)

    try
    catch error
        return VLResult(error)
    end
end

function translate_sequence(table::DataFrame, complementOperation::Function;
    logger::Union{Nothing,SimpleLogger}=nothing)


    try
    catch error
        return VLResult(error)
    end
end