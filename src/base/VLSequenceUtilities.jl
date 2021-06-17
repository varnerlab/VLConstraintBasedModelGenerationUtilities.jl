function build_transcription_reaction_table(gene_name::String, sequence::BioSequences.LongSequence; 
    polymeraseSymbol::Symbol=:RNAP, ecnumber::String="2.7.7.6", logger::Union{Nothing,SimpleLogger}=nothing)::VLResult

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
        
        # setup reactions -
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

function build_translation_reaction_table(protein_name::String, sequence::BioSequences.LongAminoAcidSeq; 
    ribosomeSymbol::Symbol=:RIBOSOME, ecnumber::String="3.1.27.10", logger::Union{Nothing,SimpleLogger}=nothing)::VLResult

    try

        # initailize -
        id_array = String[]
        forward_reaction_string = String[]
        reverse_reaction_string = String[]
        reversibility_array = Bool[]
        ec_number_array = Union{Missing,String}[]
        ribosomeSymbol = string(ribosomeSymbol)
        total_residue_count = 0
        protein_aa_map = Dict{BioSymbol, Int64}()
        aa_biosymbol_array = [
            AA_A, AA_R, AA_N, AA_D, AA_C, AA_Q, AA_E, AA_G, AA_H, AA_I, AA_L, AA_K, AA_M, AA_F, AA_P, AA_S, AA_T, AA_W, AA_Y, AA_V
        ];

        # load AA map -
        aa_metabolite_map = TOML.parsefile(joinpath(_PATH_TO_CONFIG, "AAMap.toml"))

        # build the protein_aa_map -
        for residue in aa_biosymbol_array
            protein_aa_map[residue] = count(sequence, residue)
        end

        # total residue count -
        for residue in aa_biosymbol_array
            number_per_AA = protein_aa_map[residue]
            total_residue_count = total_residue_count + number_per_AA
        end
        
        # setup reactions -
        # mRNA + RIBOSOME <=> mRNA_RIBOSOME_closed
        push!(id_array, "$(protein_name)_binding")
        push!(forward_reaction_string, "mRNA_$(protein_name)+$(ribosomeSymbol)")
        push!(reverse_reaction_string, "mRNA_$(protein_name)_$(ribosomeSymbol)_closed")
        push!(reversibility_array, true)
        push!(ec_number_array, missing)

        # mRNA_RIBOSOME_closed => mRNA_RIBOSOME_start
        push!(id_array, "mRNA_$(protein_name)_open")
        push!(forward_reaction_string, "mRNA_$(protein_name)_$(ribosomeSymbol)_closed")
        push!(reverse_reaction_string, "mRNA_$(protein_name)_$(ribosomeSymbol)_start")
        push!(reversibility_array, false)
        push!(ec_number_array, missing)

        # mRNA_RIBOSOME_translation -
        push!(id_array, "mRNA_$(protein_name)_translation")
        push!(forward_reaction_string, "mRNA_$(protein_name)_$(ribosomeSymbol)_start+$(2*total_residue_count)*M_gtp_c+$(2*total_residue_count)*M_h2o_c")
        push!(reverse_reaction_string, "mRNA_$(protein_name)+$(ribosomeSymbol)+P_$(protein_name)+$(2*total_residue_count)*M_gdp_c+$(2*total_residue_count)*M_pi_c+$(total_residue_count)*tRNA_c")
        push!(reversibility_array, false)
        push!(ec_number_array, missing)

        # charge the tRNA -
        for residue in aa_biosymbol_array
            
            # get the key symbol -
            key_value = "AA_"*(string(residue))

            # get the number and M_* of this residue -
            number_of_AA_residue = protein_aa_map[residue]
            metabolite_symbol = aa_metabolite_map[key_value]

            # build the reaction record -
            push!(id_array, "tRNA_charging_$(metabolite_symbol)_$(protein_name)")            
            push!(forward_reaction_string, "$(number_of_AA_residue)*$(metabolite_symbol)+$(number_of_AA_residue)*M_atp_c+$(number_of_AA_residue)*tRNA_c+$(number_of_AA_residue)*M_h2o_c")
            push!(reverse_reaction_string, "$(number_of_AA_residue)*$(metabolite_symbol)_tRNA_c+$(number_of_AA_residue)*M_amp_c+$(number_of_AA_residue)*M_ppi_c")
            push!(reversibility_array, false)
            push!(ec_number_array, missing)
        end

        # package into DataFrame -
        reaction_dataframe = DataFrame(id=id_array, forward=forward_reaction_string, reverse=reverse_reaction_string, reversibility=reversibility_array, ec=ec_number_array)

        # return -
        return VLResult(reaction_dataframe)
    catch error
        return VLResult(error)
    end
end

function build_transport_reaction_table()::VLResult

    try

        # initailize -
        id_array = String[]
        forward_reaction_string = String[]
        reverse_reaction_string = String[]
        reversibility_array = Bool[]
        ec_number_array = Union{Missing,String}[]
        aa_biosymbol_array = [
            AA_A, AA_R, AA_N, AA_D, AA_C, AA_Q, AA_E, AA_G, AA_H, AA_I, AA_L, AA_K, AA_M, AA_F, AA_P, AA_S, AA_T, AA_W, AA_Y, AA_V
        ];

        # load AA map -
        aa_metabolite_map = TOML.parsefile(joinpath(_PATH_TO_CONFIG, "AAMap.toml"))

        # add exchange water reactions -
        push!(id_array, "M_h2o_c_exchange")
        push!(forward_reaction_string, "M_h2o_e")
        push!(reverse_reaction_string, "M_h2o_c")
        push!(reversibility_array, true)
        push!(ec_number_array, missing)

        # add M_ppi_e exchange -
        push!(id_array, "M_ppi_c_exchange")
        push!(forward_reaction_string, "M_ppi_e")
        push!(reverse_reaction_string, "M_ppi_c")
        push!(reversibility_array, true)
        push!(ec_number_array, missing)

        # add M_amp_c exchange -
        push!(id_array, "M_amp_c_exchange")
        push!(forward_reaction_string, "M_amp_e")
        push!(reverse_reaction_string, "M_amp_c")
        push!(reversibility_array, true)
        push!(ec_number_array, missing)

        # add M_gmp_c exchange -
        push!(id_array, "M_gmp_c_exchange")
        push!(forward_reaction_string, "M_gmp_e")
        push!(reverse_reaction_string, "M_gmp_c")
        push!(reversibility_array, true)
        push!(ec_number_array, missing)

        # add M_cmp_c exchange -
        push!(id_array, "M_cmp_c_exchange")
        push!(forward_reaction_string, "M_cmp_e")
        push!(reverse_reaction_string, "M_cmp_c")
        push!(reversibility_array, true)
        push!(ec_number_array, missing)

        # add M_ump_c exchange -
        push!(id_array, "M_ump_c_exchange")
        push!(forward_reaction_string, "M_ump_e")
        push!(reverse_reaction_string, "M_ump_c")
        push!(reversibility_array, true)
        push!(ec_number_array, missing)


        # add M_atp_c exchange -
        push!(id_array, "M_atp_c_exchange")
        push!(forward_reaction_string, "M_atp_e")
        push!(reverse_reaction_string, "M_atp_c")
        push!(reversibility_array, true)
        push!(ec_number_array, missing)

        # add M_gtp_c exchange -
        push!(id_array, "M_gtp_c_exchange")
        push!(forward_reaction_string, "M_gtp_e")
        push!(reverse_reaction_string, "M_gtp_c")
        push!(reversibility_array, true)
        push!(ec_number_array, missing)

        # add M_ctp_c exchange -
        push!(id_array, "M_ctp_c_exchange")
        push!(forward_reaction_string, "M_ctp_e")
        push!(reverse_reaction_string, "M_ctp_c")
        push!(reversibility_array, true)
        push!(ec_number_array, missing)

        # add M_utp_c exchange -
        push!(id_array, "M_utp_c_exchange")
        push!(forward_reaction_string, "M_utp_e")
        push!(reverse_reaction_string, "M_utp_c")
        push!(reversibility_array, true)
        push!(ec_number_array, missing)

        
        push!(id_array, "tRNA_c_exchange")
        push!(forward_reaction_string, "tRNA_e")
        push!(reverse_reaction_string, "tRNA_c")
        push!(reversibility_array, true)
        push!(ec_number_array, missing)

        # transfer AAs -
        for residue in aa_biosymbol_array
            
            # get the key symbol -
            key_value = "AA_"*(string(residue))

            # metabilite -
            metabolite_symbol_c = aa_metabolite_map[key_value]
            metabolite_symbol_e = replace(metabolite_symbol_c, "_c"=>"_e")

            # write the record -
            push!(id_array, "$(metabolite_symbol_c)_exchange")
            push!(forward_reaction_string, metabolite_symbol_e)
            push!(reverse_reaction_string, metabolite_symbol_c)
            push!(reversibility_array, true)
            push!(ec_number_array, missing)
        end

        # package into DataFrame -
        reaction_dataframe = DataFrame(id=id_array, forward=forward_reaction_string, reverse=reverse_reaction_string, reversibility=reversibility_array, ec=ec_number_array)

        # return -
        return VLResult(reaction_dataframe)
    catch error
        return VLResult(error)
    end
end

function transcribe(sequence::BioSequences.LongSequence; complementOperation::Function=!, 
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

function transcribe(table::DataFrame, complementOperation::Function=!; 
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

function translate(sequence::BioSequences.LongAminoAcidSeq, complementOperation::Function;
    logger::Union{Nothing,SimpleLogger}=nothing)::VLResult

    try
    catch error
        return VLResult(error)
    end
end

function translate(table::DataFrame, complementOperation::Function;
    logger::Union{Nothing,SimpleLogger}=nothing)::VLResult


    try

        # get the size of the table -
        (number_of_proteins, _) = size(table)
        for protein_index in number_of_proteins
            
            # get the sequence -
            protein_seq = table[protein_index, :protein_sequence]
        end

    catch error
        return VLResult(error)
    end
end