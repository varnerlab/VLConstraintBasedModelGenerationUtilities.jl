function build_transcription_reactions_for_gene_table(table::DataFrame)::VLResult

    try

    catch error
        return VLResult(error)
    end
end

function build_translation_reactions_for_protein_table(table::DataFrame)::VLResult

    try

    catch error
        return VLResult(error)
    end
end