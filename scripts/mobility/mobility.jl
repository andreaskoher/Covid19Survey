using DrWatson
quickactivate(@__DIR__)
@show projectdir()
##
using StatsBase
using StatsPlots
using DataFrames
using Underscores
using CSV
using HTTP
using JSON
using Dates

include("mobility_utils.jl")
## ============================================================================
#                    apple
#  ============================================================================
old2new_names = Dict(
    "Capital Region of Denmark" => "capital",
    "Central Denmark Region" => "central",
    "North Denmark Region" => "north",
    "Region Zealand" => "zealand",
    "Region of Southern Denmark" => "south"
)
regions = values(old2new_names)

apple_dk = @_ projectdir("data/mobility/apple/applemobilitytrends-2022-03-29.csv") |>
    CSV.read(__, DataFrame) |>
    filter(:region => ==("Denmark") , __) |>
    DataFrames.transform(__, :region => ByRow(x -> "denmark") => :region)

apple = @_ projectdir("data/mobility/apple/applemobilitytrends-2022-03-29.csv") |>
    CSV.read(__, DataFrame) |>
    filter(:country => x->!ismissing(x) , __) |>
    # filter(:geo_type => ==("country/region"), __) |>
    filter(:country => ==("Denmark") , __) |>
    filter(:geo_type => ==("sub-region") , __) |>
    # filter(:region => x-> x ∈ keys(old2new_names) , __) |>
    DataFrames.transform(__, :region => ByRow(x -> old2new_names[x]) => :region) |>
    vcat(__, apple_dk) |>
    filter(:transportation_type => !=("transit") , __) |> #WARNING TRANSIT REMOVED
    DataFrames.select(__, Not(["geo_type", "alternative_name", "sub-region", "country"])) |>
    DataFrames.stack(__, Not([:region, :transportation_type]), variable_name="date", value_name="mobility") |>
    unstack(__, :transportation_type, :mobility) |>
    DataFrames.select(__, :date => ByRow(Date) => :date, Not(:date), AsTable([:driving,:walking]) => ByRow(x->custom_mean(x)) => :apple) |> #WARNING TRANSIT REMOVED
    rename(__, [:driving,:walking] .=> [:apple_driving,:apple_walking]) |> #WARNING TRANSIT REMOVED
    filter( Date("2020-02-15") <= _.date, __) # NOTE filter according to google


## ============================================================================
# plot
begin
    regions = unique(apple.region)
    plts = Vector{Plots.Plot}()
    for region in regions
        regional = filter(:region => ==(region), apple)
        @assert all(regional.date .== regional.date[1]:Day(1):regional.date[end])

        p = plot(title="$region", ticks=:native)
        plot!(regional.date, regional.apple_driving, lab="driving")
        plot!(regional.date, regional.apple_walking, lab="walking")
        # !all(ismissing.(regional.apple_transit)) && plot!(regional.date, regional.apple_transit, lab="transit")
        push!(plts, p)
    end
    p = plot(plts..., layout=(length(plts), 1), size=(1000, length(plts)*250), legend = true)
    p
end
## ============================================================================
#                    google
#  ============================================================================
mobility_names = [
    "retail_and_recreation_percent_change_from_baseline",
    "grocery_and_pharmacy_percent_change_from_baseline",
    "workplaces_percent_change_from_baseline",
    "residential_percent_change_from_baseline",
    "transit_stations_percent_change_from_baseline",
    # "parks_percent_change_from_baseline"
]

new_mobility_names = [
    "google_retail",
    "google_grocery",
    "google_workplaces",
    "google_residential",
    "google_transit",
    # "google_parks"
]

google_dk = @_ projectdir("data/mobility/google/Global_Mobility_Report.csv") |>
    CSV.File(__) |>
    DataFrame |>
    filter(:country_region => ==("Denmark"), __) |>
    filter( :sub_region_1 => ismissing, __) |>
    DataFrames.transform(__, :sub_region_1 => ByRow(x -> "denmark") => :region) |>
    DataFrames.transform(__, :region => disallowmissing => :region)


google = @_ projectdir("data/mobility/google/Global_Mobility_Report.csv") |>
    CSV.File(__) |>
    DataFrame |>
    filter(:country_region => ==("Denmark"), __) |>
    filter( :sub_region_1 => x->!ismissing(x), __) |>
    DataFrames.transform(__, :sub_region_1 => disallowmissing => :region) |>
    filter(:sub_region_1 => x-> x ∈ keys(old2new_names) , __) |>
    filter(:sub_region_2 => ismissing, __) |>
    DataFrames.select(__, All(), :region => ByRow(x -> old2new_names[x]) => :region) |>
    vcat(__, google_dk) |>
    DataFrames.select(__, :date, :region, mobility_names .=> new_mobility_names) |>
    DataFrames.select(__, All(), AsTable([:google_retail,:google_grocery,:google_workplaces, :google_transit]) => ByRow(mean) => :google) |>
    DataFrames.DataFrames.transform(__, [new_mobility_names; "google"] .=> ByRow(x -> x + 100); renamecols=false)

## ============================================================================
# plot
begin
    regions = unique(apple.region)
    plts = Vector{Plots.Plot}()
    for region in regions
        regional = @_ google |>
            filter(:region => ==(region), __) |>
            filter(:date => >=(Date("2020-08-01")), __) |>
            filter(:date => <=(Date("2021-02-01")), __)
        @assert all(regional.date .== regional.date[1]:Day(1):regional.date[end])

        p = plot(title="$region")
        plot!(regional.date, regional.google_retail, lab="google_retail")
        plot!(regional.date, regional.google_grocery, lab="google_grocery")
        plot!(regional.date, regional.google_workplaces, lab="google_workplaces")
        plot!(regional.date, regional.google_transit, lab="google_transit")
        # plot!(regional.date, regional.google_parks, lab="google_parks")
        # plot!(regional.date, regional.google, lab="non-residential")
        # !all(ismissing.(regional.apple_transit)) && plot!(regional.date, regional.apple_transit, lab="transit")
        push!(plts, p)
    end
    p = plot(plts..., layout=(length(plts), 1), size=(1000, length(plts)*250), legend = false, link=:x)
    p
end
## ============================================================================
#                    telco
#  ============================================================================
kom2reg = kommune2region(
      old2new_names
    ; google_fname = projectdir("data/mobility/google/Global_Mobility_Report.csv")
)

telco_raw = @_ projectdir("data/mobility/telco/telco_data.json") |>
    JSON.parsefile |>
    json2long(__ ; dateformat = dateformat"yyyy-mm-dd HH:MM:SS") |>
    DataFrames.transform(__, :source => ByRow(x->kom2reg[x]) => :source) |>
    DataFrames.transform(__, :target => ByRow(x->kom2reg[x]) => :target)

telco_dk = @_ telco_raw |>
    DataFrames.select(__, :date, :flow => :telco, :source => ByRow(x -> "denmark") => :region) |>
    groupby(__, [:date, :region]) |>
    combine(__, :telco => sum => :telco)

telco = @_ telco_raw |>
        filter(_.source == _.target, __) |>
        DataFrames.select(__, :date, :source => :region, :flow => :telco) |>
        groupby(__, [:date, :region]) |>
        combine(__, :telco => sum => :telco) |>
        vcat(__, telco_dk) |>
        sort(__, [:region, :date])
## ============================================================================
# plot
begin
    regions = unique(telco.region)
    plts = Vector{Plots.Plot}()
    for region in regions
        regional = @_ telco |>
            filter(:region => ==(region), __) |>
            filter(:date => >=(Date("2020-08-01")), __) |>
            filter(:date => <=(Date("2021-02-01")), __)
        @assert all(regional.date .== regional.date[1]:Day(1):regional.date[end])

        p = plot(title="$region")
        plot!(regional.date, regional.telco, lab="telco")
        push!(plts, p)
    end
    p = plot(plts..., layout=(length(plts), 1), size=(1000, length(plts)*250), legend = false, link=:x)
    # firefox(p; fname="apple.html")
    p
end

## ============================================================================
#                    combine google, apple and telco
#  ============================================================================
mobility_raw = @_ outerjoin(apple, google; on=[:date,:region]) |>
    outerjoin(__, telco; on=[:date,:region]) |>
    sort(__, [:region, :date]) |>
    filter( Date("2020-06-01") <= _.date <= Date("2021-03-01"), __)
# CSV.write(projectdir("data/mobility/", "mobility_raw.csv"), mobility_raw)
# mobility_raw = CSV.read(projectdir("data/mobility/", "mobility_raw.csv"), DataFrame)

## ============================================================================
#                    smooth and clean up raw data
#  ============================================================================

mobility = postprocessing_mobilty(mobility_raw)
# CSV.write(projectdir("data/mobility/", "mobility.csv"), mobility)

## plot
let
    ps = []
    for region in unique(mobility.region)
        regional = filter(:region => ==(region), mobility)
        p = plot(title="$region", xticks=:native)
        plot!(regional.date, regional.apple, lab="apple")
        plot!(regional.date, regional.google, lab="google")
        plot!(regional.date, regional.telco, lab="telco")
        push!(ps,p)
    end
    nplots = length(ps)
    plot(ps..., layout=(nplots,1), sharex=true, link=:x, size=(1000, nplots*250), legend=false, ticks=:native)
end
