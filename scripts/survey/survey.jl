using DrWatson
quickactivate(@__DIR__)
@show projectdir()
##
using CSV#
using DataFrames
using Dates
using StatsBase
using StatsPlots
using Underscores
const dateform = "yyyymmdd"
plotlyjs()
using ColorSchemes
colors = ColorSchemes.tableau_10
include( "survey_utils.jl" )

## =========================================================================
# read from original HOPE survey #NOTE full data will not be provided

# dan2eng = Dict(
#     "Hovedstaden" => "capital",
#     "Sj�lland" => "zealand",
#     "Midtjylland" => "central",
#     "Syddanmark" => "south",
#     "Nordjylland" => "north"
# )
#
# survey = @_ normpath( homedir(), "data/covidsurvey/rawcontacts.csv") |>
#     CSV.read(__, DataFrame, stringtype=String) |>
#     transform(__, :region => ByRow(x->dan2eng[x]) => :region) |>
#     select(__,
#     :region,
#     :Timings_yyyymmdd => ByRow(x -> Date("$x", dateformat"yyyymmdd") ) => :date,
#     :Q4a_1 => :family,
#     :Q4a_2 => :colleagues,
#     :Q4a_3 => :friends,
#     :Q4a_4 => :strangers,
#     :Q1_2 => :Q1_threat_to_society,
#     :Q3_2 => :Q3_avoid_contacts,
#     :Q3_4 => :Q3_avoid_vulnerable,
#     :Q3_5 => :Q3_social_distancing,
#     :Q3_6 => :Q3_avoid_crowded_places,
#     :Q3_7 => :Q3_avoid_contact_activities
#     ) |>
#     sort(__, [:region, :date])
#
# canonicaltimes(df) =  df.date[1] : Day(1) : df.date[end]
# function hascanonicaltimes(df)
#     length(canonicaltimes(df)) == length(df.date) && (return df)
#     ds = canonicaltimes(df)
#     println("missing dates: $([d for d in ds if d ∉ df.date])")
#     @assert false "missing dates"
#     df
# end
#
# df = filter(:region => ==("capital"), survey)
#
# for region in unique(survey.region)
#     @show region
#     df = filter(:region => ==(region), survey)
#     dates = unique(df.date)
#     canonical = dates[1] : Day(1) : dates[end]
#     @assert all(dates .== canonical)
# end
#
# CSV.write(projectdir("data/survey/selected_survey_responses.csv"), survey)
## =========================================================================================================
kwargs = (
    ; fname = projectdir( "data/survey/selected_survey_responses.csv")
    , outlier_threshold = (
        ; family=30 #50 #30
        , colleagues=50 #100 #50
        , friends=50 #100 #50
        , strangers=100 #1000 #100
        )
    )
survey = survey_predictors(; kwargs...)
# CSV.write(projectdir("data/survey/survey.csv"), survey)
## ============================================================================
# plot
begin
    regions = unique(survey.region)
    plts = Vector{Plots.Plot}()
    for region in regions
        regional = filter(:region => ==(region), survey)
        @assert all(regional.date .== regional.date[1]:Day(1):regional.date[end])

        p = plot(title="$region", ticks=:native)
        plot!(regional.date, regional[:,"friends-above4"], lab="friends")
        plot!(regional.date, regional[:,"strangers-above5"], lab="strangers")
        plot!(regional.date, regional[:,"colleagues-above4"], lab="colleagues")
        plot!(regional.date, regional[:,"family-above5"], lab="family")
        push!(plts, p)
    end
    p = plot(plts..., layout=(length(plts), 1), size=(1000, length(plts)*250), legend = true)
    p
end

## ==========================================================================
# quantiles and the corresponding threshold number of contacts

prelockdown = @_  rawcontacts(; kwargs...) |>
    filter(:date => >=(Date("2020-08-01")), __) |>
    filter(:date => <=(Date("2020-12-06")), __)

q = 0.85                            #0.5 #0.6 #0.7 #0.75 #0.8 #0.85 #0.9
@show quantile(prelockdown.family, q)     #0   #1   #2   #2    #3   #4    #5
@show quantile(prelockdown.friends, q)    #0   #1   #2   #3    #4   #5    #6
@show quantile(prelockdown.colleagues, q) #0   #0   #2   #3    #4   #5    #6
@show quantile(prelockdown.strangers, q)  #0   #0   #2   #3    #5   #6    #10
@show quantile(prelockdown.total, q)      #5   #7   #10  #12   #15  #19   #25

pgfplotsx()
scalefontsizes()
scalefontsizes(2)
p = let  #log10
    p1 = plot(; yformatter=:scientific, ylabel="probability", xlabel="reported contacts", xlims=(0,Inf), ylims=(0,Inf))#, title="A", titlelocation=:left, titlefontfamily="Courier Bold")
    histogram!(prelockdown.total; bins=0:220, normalize=:probability, fc=colors[1], lc=colors[1], lab=nothing, yscale = :identity )

    p2 = plot(; yformatter=:scientific, xlabel="reported contacts", xlims=(0,Inf), ylims=(1e-5,Inf), yticks = [0.1, 0.01, 0.001, 0.0001])#, title="B", titlelocation=:left, titlefontfamily="Courier Bold")
    histogram!(prelockdown.total; bins=0:220, normalize=:probability, fc=colors[1], lc=colors[1], lab=nothing, yscale = :log10 )

    plot(p1, p2, layout=(1,2), dpi=300, size=(1000, 400))
end
sum(prelockdown.total .== 0) / size(prelockdown,1) * 100
