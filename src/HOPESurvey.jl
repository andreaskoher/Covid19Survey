module Covid19Survey

# using DrWatson
# using CSV
# using MCMCChains
# using DistributionsAD
using DrWatson
using BSON
using Turing
using Distributions
using Underscores
using Random: AbstractRNG
using StatsBase
using DataFrames
using Dates
using PyCall

# postprocessing
using StatsPlots
using StatsPlots.Plots.PlotMeasures
import PyPlot as plt
using PrettyTables
using OrderedCollections
using CSV
import Pandas
# using StatsBase
# using Underscores
# using Base: @kwdef


include("utils.jl")
include("data.jl")
include("models.jl")
include("postprocessing.jl")
include("visualizations.jl")

# const regions = ["capital", "zealand", "south", "central", "north"]
end # module
