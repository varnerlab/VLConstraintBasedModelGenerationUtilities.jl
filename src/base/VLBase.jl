import Base.+

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