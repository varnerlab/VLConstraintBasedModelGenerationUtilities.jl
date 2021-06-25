using Base: has_nondefault_cmd_flags
import Base.+
import Base.!
import Base.count

function !(nucleotide::BioSymbols.DNA)::BioSymbols.RNA

    # ok: use watson-crick base pairing rules to return the opposite nucleotide
    if nucleotide == DNA_T
        return RNA_A
    elseif nucleotide == DNA_A
        return RNA_U
    elseif nucleotide == DNA_C
        return RNA_G
    elseif nucleotide == DNA_G
        return RNA_C
    else
        return RNA_N
    end
end

function +(buffer::Array{String,1}, content::String; 
    prefix::Union{String,Nothing}=nothing,suffix::Union{String,Nothing}=nothing)
    
    # create a new content line -
    new_line = content
    
    # prefix -
    if (prefix !== nothing)
        new_line = prefix*new_line
    end

    # suffix -
    if (suffix !== nothing)
        new_line = new_line*suffix
    end
    
    # cache -
    push!(buffer,new_line)
end

function read_file_from_path(path_to_file::String)::Array{String,1}

    # initialize -
    buffer = String[]

    # Read in the file -
    open("$(path_to_file)", "r") do file
        for line in eachline(file)
            +(buffer,line)
        end
    end

    # return -
    return buffer
end

function +(buffer::Array{String,1}, content_array::Array{String,1})
    for line in content_array
        push!(buffer, line)
    end
end

function check(result::VLResult)::(Union{Nothing,T} where T <: Any)

    # ok, so check, do we have an error object?
    # Yes: log the error if we have a logger, then throw the error. 
    # No: return the result.value

     # Error case -
    if (isa(result.value, Exception) == true)
        
        # get the error object -
        error_object = result.value

        # get the error message as a String -
        error_message = sprint(showerror, error_object, backtrace())
        @error(error_message)

        # throw -
        throw(result.value)
    end

    # default -
    return result.value
end

function count(sequence::BioSequences.LongSequence, bioSymbol::BioSymbol)::Int64

    # initialize -
    number_of_biosymbols = 0
    
    # iterate -
    for test_symbol in sequence
        if (test_symbol == bioSymbol)
            number_of_biosymbols = number_of_biosymbols + 1
        end
    end

    # return -
    return number_of_biosymbols
end