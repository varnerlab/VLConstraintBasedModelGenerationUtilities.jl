using Documenter
using VLConstraintBasedModelGenerationUtilities

makedocs(
    sitename = "Model",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [VLConstraintBasedModelGenerationUtilities],
    pages = [
        "Home" => "index.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/varnerlab/VLConstraintBasedModelGenerationUtilities.jl.git",
    devurl = "stable",
    devbranch = "main,
)
