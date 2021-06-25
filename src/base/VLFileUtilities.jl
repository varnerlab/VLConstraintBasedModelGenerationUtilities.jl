# === PRIVATE METHODS BELOW HERE =================================================================== #
function _build_vff_reaction_table(path_to_reaction_file::String)::DataFrame

    # vff: CSV format with record structure
    # name,[ec;...],reactant,product,reversible 

    # initialize -
    id_array = String[]
    forward_reaction_string = String[]
    reverse_reaction_string = String[]
    reversibility_array = Bool[]
    ec_number_array = Union{Missing,String}[]

    # get file buffer (array of strings) -
    vff_file_buffer = read_file_from_path(path_to_reaction_file)

    # process -
    for (_, reaction_line) in enumerate(vff_file_buffer)
        
        # skip comments and empty lines -
        if (occursin("//", reaction_line) == false && isempty(reaction_line) == false)
            
            # split around ,
            reaction_record_components_array = split(reaction_line, ",")

            # process each of the components -
            
            # 1: id -
            push!(id_array, string(reaction_record_components_array[1]))
            
            # 2: ec numbers -
            ec_number_component = reaction_record_components_array[2]
            if (ec_number_component == "[]")
                push!(ec_number_array, missing)
            else
                push!(ec_number_array, string(ec_number_component[2:end - 1]))
            end

            # 3: L phrase -
            push!(forward_reaction_string, string(reaction_record_components_array[3]))

            # 4: R phrase -
            push!(reverse_reaction_string, string(reaction_record_components_array[4]))

            # 5: reverse -
            push!(reversibility_array, parse(Bool, string(reaction_record_components_array[5])))
        end
    end

    # build the df -
    df_metabolic_reactions = DataFrame(id=id_array, forward=forward_reaction_string, reverse=reverse_reaction_string, reversibility=reversibility_array, ec=ec_number_array)

    # return -
    return df_metabolic_reactions
end

function _build_json_reaction_table(path_to_reaction_file::String)::DataFrame
end

function _build_sbml_reaction_table(path_to_reaction_file::String)::DataFrame
end

# === PRIVATE METHODS ABOVE HERE =================================================================== #

# === PUBLIC METHODS BELOW HERE ==================================================================== #
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
        sequence_data_frame = DataFrame(id=String[], description=String[], gene_sequence=Union{Missing,BioSequences.LongSequence}[])
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
        sequence_data_frame = DataFrame(id=String[], description=String[], protein_sequence=Union{Missing,BioSequences.LongSequence}[])
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

    try

        # what is the set of extensions that we can parse?
        extension_set = Set{String}()
        push!(extension_set, ".vff")
        push!(extension_set, ".json")
        push!(extension_set, ".sbml")
        push!(extension_set, ".xml")

        # get the file name from the path -
        reaction_filename = basename(path_to_reaction_file)

        # do we support this extension?
        reaction_file_extension = string(last(splitext(reaction_filename)))
        if (in(reaction_file_extension, extension_set) == false)
            throw(ArgumentError("Reaction file extension not supported: $(reaction_file_extension)"))
        end

        # ok: if we get here then we support this extension -
        reaction_table = nothing
        if (reaction_file_extension == ".vff")
            reaction_table = _build_vff_reaction_table(path_to_reaction_file)
        elseif (reaction_file_extension == ".json")
            reaction_table = _build_json_reaction_table(path_to_reaction_file)
        elseif (reaction_file_extension == ".sbml" || reaction_file_extension == ".xml")
            reaction_table = _build_sbml_reaction_table(path_to_reaction_file)
        end

        # return -
        return VLResult(reaction_table)
    catch error
        return VLResult(error)
    end
end
# === PUBLIC METHODS ABOVE HERE ==================================================================== #