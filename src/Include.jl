# setup project paths -
const _PATH_TO_SRC = dirname(pathof(@__MODULE__))
const _PATH_TO_BASE = joinpath(_PATH_TO_SRC, "base")

# packages that I'm going to use ...
using Test
using DataFrames
using CSV
using BioSequences
using FASTX
using DelimitedFiles
using Logging

# load my codes ...
include(joinpath(_PATH_TO_BASE, "VLTypes.jl"))
include(joinpath(_PATH_TO_BASE, "VLBase.jl"))
include(joinpath(_PATH_TO_BASE, "VLSequenceUtilities.jl"))
include(joinpath(_PATH_TO_BASE, "VLMetabolicUtilities.jl"))
include(joinpath(_PATH_TO_BASE, "VLFileUtilities.jl"))
