using DrWatson
quickactivate(@__DIR__)
@show projectdir()
#-----------------------------------------------------------------------------
# read input
using ArgParse
using Dates
argtable = ArgParseSettings(
    description="This script samples from Regional.model using NUTS."
)
@add_arg_table! argtable begin
    "--samples", "-s"
        help = "number of samples after warm-up"
        arg_type = Int
        default = 1000
    "--chains", "-c"
        help = "number of chains sampled with multi-threading"
        arg_type = Int
        default = 4
    "--warmup", "-w"
        help = "number of samples to use for warmup/adaptation"
        arg_type = Int
        default = 1000
    "--predictors", "-p"
        help = "choose from the predictors in data/data_with_predictors.csv. separate multiple predictors with a comma, e.g: total_above18,google,apple,telco"
        arg_type = String
        default = nothing
    "--foldername", "-f"
        arg_type = String
        default = "$(today())"
    "--thinning", "-t"
        help = "thinning"
        arg_type = Int
        default = 1
    "--tree", "-e"
        help = "max tree size"
        arg_type = Int
        default = 5
    "--model", "-m"
        arg_type = String
        default = "nonparametricmodel"
    # "--prefix", "-x"
    #     arg_type = String
    #     default = ""
end
parsed_args = parse_args(ARGS, argtable)

#-----------------------------------------------------------------------------
@info "load packages"
using Covid19Survey
using CSV
using DataFrames
using BSON
using InlineStrings
using WeakRefStrings
using PooledArrays
using Turing
using Memoization
using ReverseDiff
using StatsPlots
using PrettyTables
plotlyjs()
Turing.setrdcache(true)
setadbackend(:reversediff)
@info "number of threads available: $(Base.Threads.nthreads())"
# using Random, Dates, Turing
# using OrderedCollections
## ==========================================================================
# load data
data = CSV.read(projectdir("data/data_without_predictors.csv"), DataFrame)
kwargs = BSON.load(projectdir("data/data_without_predictors.bson"))

## ==========================================================================
@info "initialize model"

tdata = Covid19Survey.turingformat(data; kwargs...)
model = try
    getfield(Covid19Survey, Symbol(parsed_args["model"]))
catch e
    @error "it seems that your model $(parsed_args["model"]) has not been defined in the Covid19Survey module."
    throw(e)
end
m = model(tdata.turing_data)
m() #test

## ==========================================================================
@info "start sampling"

@time chain = let
    @unpack thinning, warmup, samples, chains, tree = parsed_args
    if chains > 1
        sample(m, NUTS(warmup, 0.99; max_depth=tree), MCMCThreads(), samples, chains; progress=true, thinning, save_state=true)
    else
        sample(m, NUTS(warmup, 0.99; max_depth=tree), samples; progress=true, thinning, save_state=true)
    end
end

## ==========================================================================
@info "save"
fdir = projectdir("reports", parsed_args["foldername"])
fname = Covid19Survey.save_results(fdir, parsed_args, tdata, chain)
## ==========================================================================
@info "post-processing"
Covid19Survey.postprocessing(fname)
