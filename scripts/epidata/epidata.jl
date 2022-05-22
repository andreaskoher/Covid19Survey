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

include(projectdir("scripts/epidata/epidata_utils.jl"))
## ==========================================================

kwargs = Dict(
    :conditions   => [:date => x -> x <= Date("2021-12-09")], #:cumulative_hospit => x -> x > 30
    :fname_cases  => projectdir( "data/epidata/08_bekraeftede_tilfaelde_pr_dag_pr_regions.csv" ),
    :fname_hospit => projectdir( "data/epidata/06_nye_indlaeggelser_pr_region_pr_dag.csv" ),
    :fname_deaths => projectdir( "data/epidata/07_antal_doede_pr_dag_pr_region.csv" ),
    :regions      => ["capital", "zealand", "south", "central", "north"],
    :populations  => Dict(
        "capital" => 1_855_084,
        "zealand" => 838_840,
        "south"   => 1_223_634,
        "central" => 1_332_048,
        "north"   => 590_439), # from statistics denmark
    :dan2eng      => Dict(
        "Dato"          => "date",
        "Total"         => "total",
        "Hovedstaden"   => "capital",
        "Sj\xe6lland"   => "zealand",
        "Sjælland"      => "zealand",
        "Syddanmark"    => "south",
        "Midtjylland"   => "central",
        "Nordjylland"   => "north",
        "Ukendt Region" => "unknown"
        ),
)

data = epidata(; kwargs...)

CSV.write(projectdir("data/epidata/epidata.csv"), data)

## ==========================================================
# plot

plotlyjs()
begin
    regions = unique(data.region)
    for region in regions
        regional = filter(:region => ==(region), data)

        plts = Vector{Plots.Plot}()

        p = plot(title="$region: observed cases")
        plot!(regional.date, regional.cases)
        push!(plts, p)

        p = plot(title="incidence")
        plot!(regional.date, regional.cases_incidence)
        push!(plts, p)

        p = plot(title="cases per 100.000")
        plot!(regional.date, regional.cases_per_100000)
        push!(plts, p)

        p = plot(title="hospitalizations")
        plot!(regional.date, regional.hospit)
        push!(plts, p)

        p = plot(plts..., layout=(length(plts), 1), size=(1000, 250*length(plts)), legend = false, ticks=:native)
        display(p)
    end
end
