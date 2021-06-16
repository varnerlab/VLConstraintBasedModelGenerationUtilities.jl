function build_transcription_reaction_table(gene_name::String, sequence::BioSequences.LongSequence; 
    polymeraseSymbol::Symbol=:RNAP, ecnumber::String="EC 2.7.7.6", logger::Union{Nothing,SimpleLogger}=nothing)::VLResult

    try

        # initialize -
        id_array = String[]
        forward_reaction_string = String[]
        reverse_reaction_string = String[]
        reversibility_array = Bool[]
        ec_number_array = Union{Missing,String}[]
        polymeraseSymbol = string(polymeraseSymbol)

        # get a count on G,C,A,U -    
        number_of_A = count(sequence, DNA_A)
        number_of_U = count(sequence, DNA_T)
        number_of_C = count(sequence, DNA_C)
        number_of_G = count(sequence, DNA_G)
        total_nucleotides = number_of_A + number_of_U + number_of_C + number_of_G
        
        # setup reaction -
        # gene + RNAP <=> gene_RNAP_closed
        push!(id_array, "$(gene_name)_binding")
        push!(forward_reaction_string, "$(gene_name)+$(polymeraseSymbol)")
        push!(reverse_reaction_string, "$(gene_name)_$(polymeraseSymbol)_closed")
        push!(reversibility_array, true)
        push!(ec_number_array, missing)

        # gene_RNAP_closed => gene_RNAP_open
        push!(id_array, "$(gene_name)_open")
        push!(forward_reaction_string, "$(gene_name)_$(polymeraseSymbol)_closed")
        push!(reverse_reaction_string, "$(gene_name)_$(polymeraseSymbol)_open")
        push!(reversibility_array, false)
        push!(ec_number_array, missing)

        # transcription -
        push!(id_array, "$(gene_name)_transcription")
        push!(forward_reaction_string, "$(gene_name)_$(polymeraseSymbol)_open+$(number_of_A)*M_atp_c+$(number_of_U)*M_utp_c+$(number_of_C)*M_ctp_c+$(number_of_G)*M_gtp_c+$(total_nucleotides)*M_h2o_c")
        push!(reverse_reaction_string, "mRNA_$(gene_name)+$(gene_name)+$(polymeraseSymbol)+$(total_nucleotides)*M_ppi_c")
        push!(reversibility_array, false)
        push!(ec_number_array, ecnumber)

        # mRNA degradation -
        push!(id_array, "mRNA_$(gene_name)_degradation")
        push!(forward_reaction_string, "mRNA_$(gene_name)")
        push!(reverse_reaction_string, "$(number_of_A)*M_amp_c+$(number_of_U)*M_ump_c+$(number_of_C)*M_cmp_c+$(number_of_G)*M_gmp_c")
        push!(reversibility_array, false)
        push!(ec_number_array, missing)

        # package into DataFrame -
        reaction_dataframe = DataFrame(id=id_array, forward=forward_reaction_string, reverse=reverse_reaction_string, reversibility=reversibility_array, ec=ec_number_array)

        # return -
        return VLResult(reaction_dataframe)
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

function transcribe_sequence(sequence::BioSequences.LongSequence; complementOperation::Function=!, 
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
        new_seq_from_array = LongSequence{RNAAlphabet{4}}(new_sequence_array)

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