


# === PUBLIC METHODS BELOW HERE =================================================================== #
function build_gene_table(path_to_gene_file::String; 
    logger::Union{Nothing,SimpleLogger}=nothing)::VLResult

    try

        # initialize -
        local_data_row = Union{String,Symbol,Missing,FASTA.Record}[]

        # ok - lets load the gene file
        open(FASTA.Reader, path_to_gene_file) do reader
            for record in reader
                push!(local_data_row, record)
            end
        end

        # create a data frame -
        sequence_data_frame = DataFrame(id=String[], description=String[], sequence=Union{String,Missing,BioSequences.LongSequence}[])
        sequence_data_record = Union{String,Missing,BioSequences.LongSequence}[]
        for record in local_data_row
            
            # get attributes from record -
            id_value = FASTX.FASTA.identifier(record)
            description_value = FASTX.FASTA.description(record)
            sequence_value = FASTX.FASTA.sequence(record)
            
            # pack -
            push!(sequence_data_record, id_value)
            push!(sequence_data_record, description_value)
            push!(sequence_data_record, sequence_value)
            push!(sequence_data_frame, tuple(sequence_data_record...))
        end

        # return -
        return VLResult(sequence_data_frame)
    catch error
        return VLResult(error)
    end
end

function build_gene_table(gene_path_array::Array{String,1}; 
    logger::Union{Nothing,SimpleLogger}=nothing)

    try

        # initialize -
        gene_table_dictionary::Dict{String,DataFrame}()

         # ok, so let process the array of files -
        for file_path in gene_path_array
            sequence_data_frame = build_gene_table(file_path; logger=logger) |> check
            gene_table_dictionary[file_path] = sequence_data_frame
        end

        # return a dictionary of gene tables -
        return VLResult(gene_table_dictionary)
    catch error
        return VLResult(error)
    end
end

function build_protein_table(path_to_protein_file::String; 
    logger::Union{Nothing,SimpleLogger}=nothing)::VLResult

    try

        # initialize -
        local_data_row = Union{String,Symbol,Missing,FASTA.Record}[]

        # ok - lets load the gene file
        open(FASTA.Reader, path_to_protein_file) do reader
            for record in reader
                push!(local_data_row, record)
            end
        end

        # create a data frame -
        sequence_data_frame = DataFrame(id=String[], description=String[], sequence=Union{String,Missing,BioSequences.LongSequence}[])
        for record in local_data_row

            # init a row -
            sequence_data_record = Union{String,Missing,BioSequences.LongSequence}[]
            
            # get attributes from record -
            id_value = FASTX.FASTA.identifier(record)
            description_value = FASTX.FASTA.description(record)
            sequence_value = FASTX.FASTA.sequence(record)
            
            # pack -
            push!(sequence_data_record, id_value)
            push!(sequence_data_record, description_value)
            push!(sequence_data_record, sequence_value)
            push!(sequence_data_frame, tuple(sequence_data_record...))
        end

        # return -
        return VLResult(sequence_data_frame)
    catch error
        return VLResult(error)
    end
end

function build_protein_table(protein_path_array::Array{String,1}; 
    logger::Union{Nothing,SimpleLogger}=nothing)::VLResult

    try

        # initialize -
        protein_table_dictionary::Dict{String,DataFrame}()

        # ok, so let process the array of files -
        for file_path in protein_path_array
            sequence_data_frame = build_protein_table(file_path; logger=logger) |> check
            protein_table_dictionary[file_path] = sequence_data_frame
        end

        # return a dictionary of protein tables -
        return VLResult(protein_table_dictionary)
    catch error
        return VLResult(error)
    end
end

function build_metabolic_reaction_table(path_to_reaction_file::String; 
    logger::Union{Nothing,SimpleLogger}=nothing)::VLResult

    # what is the set of extensions that we can parse?
    


    try
    catch error
        return VLResult(error)
    end
end
# === PUBLIC METHODS ABOVE HERE =================================================================== #