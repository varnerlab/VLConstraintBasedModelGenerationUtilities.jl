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
    lower_bound_array = Union{Missing,Float64}[]
    upper_bound_array = Union{Missing,Float64}[]

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
            is_reversible = parse(Bool, string(reaction_record_components_array[5]))
            push!(reversibility_array, is_reversible)

            # lower bound and upper bound -
            if (is_reversible == true)
                push!(lower_bound_array, -Inf)
            else
                push!(lower_bound_array, 0.0)
            end
            push!(upper_bound_array, Inf)
        end
    end

    # build the df -
    df_metabolic_reactions = DataFrame(id=id_array, forward=forward_reaction_string, reverse=reverse_reaction_string, reversibility=reversibility_array, 
        ec=ec_number_array, lower_bound=lower_bound_array, upper_bound=upper_bound_array)

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
    logger::Union{Nothing,SimpleLogger}=nothing)::Some

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
        return Some(sequence_data_frame)
    catch error
        
        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
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
        return Some(gene_table_dictionary)
    catch error
        
        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end

function build_protein_table(path_to_protein_file::String; 
    logger::Union{Nothing,SimpleLogger}=nothing)::Some

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
        return Some(sequence_data_frame)
    catch error
        
        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end

function build_protein_table(protein_path_array::Array{String,1}; 
    logger::Union{Nothing,SimpleLogger}=nothing)::Some

    try

        # initialize -
        protein_table_dictionary::Dict{String,DataFrame}()

        # ok, so let process the array of files -
        for file_path in protein_path_array
            sequence_data_frame = build_protein_table(file_path; logger=logger) |> check
            protein_table_dictionary[file_path] = sequence_data_frame
        end

        # return a dictionary of protein tables -
        return Some(protein_table_dictionary)
    catch error
        
        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end

function build_metabolic_reaction_table(path_to_reaction_file::String; 
    logger::Union{Nothing,SimpleLogger}=nothing)::Some

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
        return Some(reaction_table)
    catch error
        
        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end

function write_system_model_file(;path::String, stoichiometric_matrix::Array{Float64,2}, species_bounds_array::Array{Float64,2}, 
    flux_bounds_array::Array{Float64,2}, reaction_table::Union{DataFrame,Nothing}=nothing)::Some

    try

        # initialize -
        system_dictionary = Dict{String,Union{Array{Float64,2},DataFrame,Nothing}}()
        
        # pack stuff into the dictionary -
        system_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
        system_dictionary["species_bounds_array"] = species_bounds_array
        system_dictionary["flux_bounds_array"] = flux_bounds_array
        system_dictionary["reaction_table"] = reaction_table

        # write -
        bson(path, system=system_dictionary)

        # return -
        return Some(nothing)
    catch error

        # get the original error message -
        error_message = sprint(showerror, error, catch_backtrace())
        vl_error_obj = ErrorException(error_message)

        # Package the error -
        return Some(vl_error_obj)
    end
end
# === PUBLIC METHODS ABOVE HERE ==================================================================== #