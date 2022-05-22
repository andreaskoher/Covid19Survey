using DrWatson
quickactivate(@__DIR__)
@show projectdir()
##
using CSV
using DataFrames
using Dates
using StatsPlots
using BSON
using StatsBase
using Underscores
using ColorSchemes
using PrettyTables
colors = ColorSchemes.tableau_10.colors
include("data_utils.jl")
## ===========================================================================
# make dataframe
regions = ["capital", "zealand", "south", "central", "north"]
kwargs = (
	; regions
	, num_impute  = 40
    , num_predictions = 0
    , epidata_kwargs = (
        fname         = projectdir("data/epidata/epidata.csv"),
        condition    = :date => x -> Date("2020-08-01") <= x <= Date("2021-02-21")
        )

    , predictor_kwargs = (
        survey_kwargs     = (
			fname         = projectdir("data/survey/survey.csv"),
            shift         = -1),
        mobil_kwargs  = (
            fname        = projectdir("data/mobility/mobility.csv"),
            ),
        relativechange = true,
        # standardize   = false,
        predictors    = nothing,
		condition     = :date => x -> Date("2020-12-01") <= x <= Date("2021-02-01")
        )
    , rw_kwargs = (
        step      = 1,
        constant_prediction = true,
        enddate   = nothing,
        semiparam = true
        )
)


data = load_data(; kwargs...)
CSV.write(projectdir("data/data_without_predictors.csv"), data)
BSON.bson(projectdir("data/data_without_predictors.bson"); kwargs...)

## ===========================================================================
# plot data
plotlyjs()
begin
    regions = unique(data.region)
    for region in regions
        regional = filter(:region => ==(region), data)
	    bs = inference_boundaries(regional)

	    plts = Vector{Plots.Plot}()

	    p = plot(title="$region: observed cases (incidence)")
	    hover = ["$d: $i" for (d,i) in zip(string.(regional.date), regional.cases_incidence)]
	    plot!(regional.date, regional.cases_incidence, hover=hover, lw = 3, label="7-day incidence")
	    push!(plts, p)

	    p = plot(title="hospitalizations")
	    plot!(regional.date, regional.hospit, lw = 3, label="hospitalizations")
	    push!(plts, p)

		p = plot(title="daily random walk")
	    plot!(regional.date, regional.rw, lab=nothing, lw = 3, ylabel="number of steps")
	    push!(plts, p)

		nps = length(plts)
	    p = plot(plts..., layout=(nps, 1), size=(1000, 300*nps), legend = true)
		display(p)
	end
end
