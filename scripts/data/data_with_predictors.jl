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
        fname        = projectdir("data/epidata/epidata.csv"),
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
        predictors    = Not(["date", "nrow"]),
		condition     = :date => x -> Date("2020-12-01") <= x <= Date("2021-02-01")
        )
    , rw_kwargs = (
        step      = 7,
        constant_prediction = true,
        enddate   = nothing,
        semiparam = true
        )
)


data = load_data(; kwargs...)
CSV.write(projectdir("data/data_with_predictors.csv"), data)
BSON.bson(projectdir("data/data_with_predictors.bson"); kwargs...)

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
		vspan!(bs, label="parametrized Rt", α=0.2, fc=:midnightblue)
	    plot!(regional.date, regional.cases_incidence, hover=hover, lw = 3, label="7-day incidence")
	    push!(plts, p)

	    p = plot(title="hospitalizations")
		vspan!(bs, label="parametrized Rt", α=0.2, fc=:midnightblue)
	    plot!(regional.date, regional.hospit, lw = 3, label="hospitalizations")
	    push!(plts, p)

		p = plot(title="risk-taking behaviour")
		vspan!(bs, label="parametrized Rt", α=0.2, fc=:midnightblue)
		plot!(regional.date, regional[!,"total-above18"], lab="total contacts", lw = 3)
		push!(plts, p)

	    p = plot(title="context-dependent risk-taking behaviour")
		vspan!(bs, label="parametrized Rt", α=0.2, fc=:midnightblue)
		plot!(regional.date, regional[:,"strangers-above5"], lab="strangers", lw = 3)
		plot!(regional.date, regional[:,"friends-above4"], label="friends", lw = 3)
		plot!(regional.date, regional[:,"colleagues-above4"], label="colleagues", lw = 3)
		push!(plts, p)

		p = plot(title="context-dependent risk-taking behaviour (family)")
		vspan!(bs, label="parametrized Rt", α=0.2, fc=:midnightblue)
	    plot!(regional.date, regional[:,"family-above3"], lab="family", lw = 3)
		push!(plts, p)

		p = plot(title="mobility")
		vspan!(bs, label="parametrized Rt", α=0.2, fc=:midnightblue)
	    plot!(regional.date, regional.google, lab="google", lw = 3)
		plot!(regional.date, regional.apple, label="apple", lw = 3)
		plot!(regional.date, regional.telco, label="telco", lw = 3)
		push!(plts, p)

		p = plot(title="weekly random walk")
		vspan!(bs, label="parametrized Rt", α=0.2, fc=:midnightblue)
	    plot!(regional.date, regional.rw, lab=nothing, lw = 3, ylabel="number of steps")
	    push!(plts, p)

		nps = length(plts)
	    p = plot(plts..., layout=(nps, 1), size=(1000, 300*nps), legend = true)
		display(p)
	end
end
