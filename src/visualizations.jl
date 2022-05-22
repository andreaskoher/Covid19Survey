abstract type PlottingRecipe end

function firefox(p; fname="tmp.html", fdir=normpath(homedir(), ".tmp"))
    # mkpath(dir)
    filepath = normpath(fdir, fname)
    savefig(p, filepath)
    run(`firefox $(filepath)`, wait=false)
    return p
end

## ===========================================================================
#              plot timeseries with confidence band
# ============================================================================

"""
    plot_confidence_timeseries!(p::Plot, data::AbstractMatrix{<:Real}; label="", kwargs...)

Plots confidence intervals for the time-series represented by `data`.
Assumes each row corresponds to the samples for a single time-step.
"""
function plot_confidence_timeseries!(p::Plots.AbstractPlot, data::AbstractMatrix; no_label = false, label="", color=:peru, α=1., lw=2, kwargs...)
    intervals = [0.025, 0.25, 0.5, 0.75, 0.975]

    qs = [quantile(v, intervals) for v in eachrow(data)]
    m  = [mean(v) for v in eachrow(data)]
    llq, lq, mq, uq, uuq = (eachrow(hcat(qs...))..., )
    plot!(m; ribbon=(mq - llq, uuq - mq), color, α, lw, label=(no_label ? "" : "$(label) (95% quantiles)"), hover=mq, kwargs...)

    return p
end

"""
    plot_confidence_timeseries!(p::Plot, time::AbstractVector, data::AbstractMatrix{<:Real}; label="", kwargs...)

Plots confidence intervals for the time-series represented by `data` and a vector of Dates.
Assumes each row corresponds to the samples for a single time-step.
"""
function plot_confidence_timeseries!(p::Plots.AbstractPlot, time::AbstractVector, data::AbstractMatrix; no_label = false, label="", α=1., lw=2, color=:peru, kwargs...)
    intervals = [0.025, 0.25, 0.5, 0.75, 0.975]

    qs = [quantile(v, intervals) for v in eachrow(data)]
    m  = [mean(v) for v in eachrow(data)]
    llq, lq, mq, uq, uuq = (eachrow(hcat(qs...))..., )
    strtime = time .|> Symbol .|> String
    #plot!(time, mq, ribbon=(mq - llq, uuq - mq), c=c, α=0.5, linewidth=0, label=(no_label ? "" : "$(label) (95% quantiles)"), hover=strtime, kwargs...)
    plot!(time, m; ribbon=(mq - llq, uuq - mq), color, α, lw, label=(no_label ? "" : "$(label) (95% quantiles)"), hover=mq, kwargs...)

    return p
end

"""
    plot_confidence_timeseries(data::AbstractMatrix{<:Real}; label="", kwargs...)

See `plot_confidence_timeseries!`.
"""
plot_confidence_timeseries(data::AbstractVector; kwargs...) = plot_confidence_timeseries!(plot(), hcat(data...) ; kwargs...)
plot_confidence_timeseries(data::AbstractMatrix; kwargs...) = plot_confidence_timeseries!(plot(), data; kwargs...)
plot_confidence_timeseries!(p::Plots.AbstractPlot, data::AbstractVector; kwargs...) = plot_confidence_timeseries!(p, hcat(data...) ; kwargs...)

"""
    plot_confidence_timeseries(data::AbstractMatrix{<:Real}; label="", kwargs...)

See `plot_confidence_timeseries!`.
"""
plot_confidence_timeseries(
    time::AbstractVector,
    data::AbstractVector;
    kwargs...
) = plot_confidence_timeseries!(plot(), time, hcat(data...); kwargs...)
plot_confidence_timeseries!(
    p   ::Plots.AbstractPlot,
    time::AbstractVector,
    data::AbstractVector;
    kwargs...
) = plot_confidence_timeseries!(p, time, hcat(data...); kwargs...)
plot_confidence_timeseries(
      time::AbstractVector
    , data::AbstractMatrix; kwargs...
) = plot_confidence_timeseries!(plot(), time, data; kwargs...)

"""
    plot_confidence_timeseries!(p::Plot, data::AbstractMatrix{<:Real}; label="", kwargs...)

Plots confidence intervals for the time-series represented by `data`.
Assumes each row corresponds to the samples for a single time-step.
"""
function plot_confidence_timeseries!(p::Plots.AbstractPlot, data::AbstractMatrix; no_label = false, label="", color=:peru, α=1., lw=2, kwargs...)
    intervals = [0.025, 0.25, 0.5, 0.75, 0.975]

    qs = [quantile(v, intervals) for v in eachrow(data)]
    m  = [mean(v) for v in eachrow(data)]
    llq, lq, mq, uq, uuq = (eachrow(hcat(qs...))..., )
    plot!(m; ribbon=(mq - llq, uuq - mq), color, α, lw, label=(no_label ? "" : "$(label) (95% quantiles)"), hover=mq, kwargs...)

    return p
end

"""
    plot_confidence_timeseries!(p::Plot, time::AbstractVector, data::AbstractMatrix{<:Real}; label="", kwargs...)

Plots confidence intervals for the time-series represented by `data` and a vector of Dates.
Assumes each row corresponds to the samples for a single time-step.
"""
function plot_confidence_timeseries!(p::Plots.AbstractPlot, time::AbstractVector, data::AbstractMatrix; no_label = false, label="", α=1., lw=2, color=:peru, kwargs...)
    intervals = [0.025, 0.25, 0.5, 0.75, 0.975]

    qs = [quantile(v, intervals) for v in eachrow(data)]
    m  = [mean(v) for v in eachrow(data)]
    llq, lq, mq, uq, uuq = (eachrow(hcat(qs...))..., )
    strtime = time .|> Symbol .|> String
    #plot!(time, mq, ribbon=(mq - llq, uuq - mq), c=c, α=0.5, linewidth=0, label=(no_label ? "" : "$(label) (95% quantiles)"), hover=strtime, kwargs...)
    plot!(time, m; ribbon=(mq - llq, uuq - mq), color, α, lw, label=(no_label ? "" : "$(label) (95% quantiles)"), hover=mq, kwargs...)

    return p
end

"""
    plot_confidence_timeseries(data::AbstractMatrix{<:Real}; label="", kwargs...)

See `plot_confidence_timeseries!`.
"""
plot_confidence_timeseries(data::AbstractVector; kwargs...) = plot_confidence_timeseries!(plot(), hcat(data...) ; kwargs...)
plot_confidence_timeseries(data::AbstractMatrix; kwargs...) = plot_confidence_timeseries!(plot(), data; kwargs...)
plot_confidence_timeseries!(p::Plots.AbstractPlot, data::AbstractVector; kwargs...) = plot_confidence_timeseries!(p, hcat(data...) ; kwargs...)

"""
    plot_confidence_timeseries(data::AbstractMatrix{<:Real}; label="", kwargs...)

See `plot_confidence_timeseries!`.
"""
plot_confidence_timeseries(
    time::AbstractVector,
    data::AbstractVector;
    kwargs...
) = plot_confidence_timeseries!(plot(), time, hcat(data...); kwargs...)
plot_confidence_timeseries!(
    p   ::Plots.AbstractPlot,
    time::AbstractVector,
    data::AbstractVector;
    kwargs...
) = plot_confidence_timeseries!(p, time, hcat(data...); kwargs...)
plot_confidence_timeseries(
      time::AbstractVector
    , data::AbstractMatrix; kwargs...
) = plot_confidence_timeseries!(plot(), time, data; kwargs...)

Plots.plot(pr::PlottingRecipe, args...; kwargs...) = plot!(plot(), pr, args...; kwargs...)


# ============================================================================
#              plotting recipees
# ============================================================================

function find_boundaries(xs::AbstractVector)
    changepoints = diff(xs)
    first(xs) == false ? pushfirst!(changepoints, 0) : pushfirst!(changepoints, 1)
    starts = findall(==(1), changepoints)
    ends = findall(==(-1), changepoints) .- 1
    length(starts) > length(ends) && push!(ends, length(xs))
    return (;starts, ends)
end

function inference_boundaries(df::DataFrame)
    bounds = find_boundaries(df.inference)
    startdates = df.date[bounds.starts]
    enddates = df.date[bounds.ends]
    bounds = [[s,e] for (s,e) in zip(startdates, enddates)]
    bounds = vcat(bounds...)
    return bounds
end


# =============================================================================
# for a single region

struct ObservationPlottingRecipe{Tsd, Ted, To, Te, Tp, Tl} <: PlottingRecipe
    startdate ::Tsd
    enddate   ::Ted
    observed  ::To
    expected  ::Te
    predicted ::Tp
    label     ::Tl
end

struct ObservationWithPredictorsPlottingRecipe{Tsd, Ted, To, Te, Tp, Tl, Psd} <: PlottingRecipe
    startdate ::Tsd
    enddate   ::Ted
    observed  ::To
    expected  ::Te
    predicted ::Tp
    label     ::Tl
    pstartdate::Psd
end

function expected(data, gp, region, label)
    @unpack dates, num_impute, regions = data
    # @assert label in ["cases", "hospitalizations", "deaths", "infections"]
    i  = region isa Integer ? region : findfirst(==(region), regions)
    dates = data.dates[i][num_impute+1:end] #data.turing_data.dates_turing[i]
    if label == "expected"
        values = gp.expecteds[i]
        return (; dates, values )
    elseif label == "newly infected"
        n = length(gp.newly_infecteds[i][1])
        values = getindex.(gp.newly_infecteds[i], Ref(num_impute+1:n))
        return (; dates, values )
    else
        @error "unkown label $label. Choose from [expected, newly infected]"
        throw(KeyError)
    end
end

function observed(data, region)
    i = region isa Integer ? region : findfirst(==(region), data.regions)
    r = region isa Integer ? regions[region] : region
    regional = data.data[r]
    dates = regional.date
    values = regional[:, data.observables]
    return (; dates, values )
end

function predicted(data, gp, region)
    @unpack dates, num_impute, regions = data
    i  = region isa Integer ? region : findfirst(==(region), regions)
    dates = data.dates[i][num_impute+1:end]
    values = gp.predicted[i]
    return (; dates, values )
end

function startdate(data::NamedTuple, region)
    # @assert label in ["cases", "hospitalizations", "deaths"]
    i = region isa Integer ? region : findfirst(==(region), data.regions)
    return data.observationsstart[i]
    # s = if label == "cases"
    #     data.turing_data.casemodel.starts[i]
    # elseif label == "hospitalizations"
    #     data.turing_data.hospitmodel.starts[i]
    # else
    #     data.turing_data.deathmodel.starts[i]
    # end
    # return data.dates[i][s]
end

function enddate(data::NamedTuple, region)
    # @assert label in ["cases", "hospitalizations", "deaths"]
    i = region isa Integer ? region : findfirst(==(region), data.regions)
    return data.observationsend[i]
    # s = if label == "cases"
    #     data.turing_data.casemodel.stops[i]
    # elseif label == "hospitalizations"
    #     data.turing_data.hospitmodel.stops[i]
    # else
    #     data.turing_data.deathmodel.stops[i]
    # end
    # return data.dates[i][s]
end

function ObservationPlottingRecipe(data::NamedTuple, gp, region, label)
    e = expected(data, gp, region, label)
    o = observed(data, region)
    p = predicted(data, gp, region)
    sd = startdate(data, region)
    ed = enddate(data, region)
    isnothing(data.predictors) && return ObservationPlottingRecipe( sd, ed, o, e, p, label)

    pstartdate = inference_boundaries(data.data[region])
    return ObservationWithPredictorsPlottingRecipe( sd, ed, o, e, p, label, pstartdate)
end

function Plots.plot!(p::Plots.Plot, r::ObservationPlottingRecipe; plot_only_fit=false, ylabel = "counts", markersize=5,
    color_predicted = :midnightblue, color_expected = :midnightblue, color_observed = :black,
    kwargs...
)
    o  = r.observed
    e  = r.expected
    pr = r.predicted
    ed = r.enddate
    sd = r.startdate

    if !plot_only_fit
        !isnothing(ed) && vline!(p, [ed], lab="end observations", lw=2, lc=:black, hover="$ed")
        !isnothing(ed) && vline!(p, [sd], lab="start observations", lw=2, lc=:black, hover="$sd", ls=:dash)
    end
    plot_confidence_timeseries!(p, pr.dates, pr.values; label = "predicted $(r.label)", color=color_predicted, kwargs...) #Dict(hover=>strdates)
    plot_confidence_timeseries!(p, e.dates, e.values; label = "expected $(r.label)", color=color_expected, kwargs...) #Dict(hover=>strdates)
    scatter!(p, o.dates, o.values; α=1., lc=:match, lab="observed $(r.label)", color=color_observed, markersize, ylabel, msw=0.)
    if plot_only_fit
        xlo = Dates.value.(sd) .- 1
        xup = Dates.value.(ed) .+ 1
        yup = maximum( o.values[sd .< o.dates .< ed] ) * 1.1
        xlims!(p, xlo, xup)
        ylims!(p, 0, yup)
    end
    return p
end

function Plots.plot!(p::Plots.Plot, r::ObservationWithPredictorsPlottingRecipe; plot_only_fit=false, ylabel = "counts", markersize=5,
    color_predicted = :lightsteelblue, color_expected = :midnightblue, color_observed = :black,
    kwargs...
)
    o  = r.observed
    e  = r.expected
    pr = r.predicted
    ed = r.enddate
    sd = r.startdate
    pd = r.pstartdate

    if !plot_only_fit
        !isnothing(ed) && vline!(p, [ed], lab="end observations", lw=2, lc=:black, hover="$ed")
        !isnothing(ed) && vline!(p, [sd], lab="start observations", lw=2, lc=:black, hover="$sd", ls=:dash)
        vspan!(p, pd; label="fit to survey data", α=0.2, fc=:lightgray)
    end
    plot_confidence_timeseries!(p, pr.dates, pr.values; label = "predicted $(r.label)", color=color_predicted, kwargs...) #Dict(hover=>strdates)
    !isnothing(color_expected) && plot_confidence_timeseries!(p, e.dates, e.values; label = "expected $(r.label)", color=color_expected, kwargs...) #Dict(hover=>strdates)
    scatter!(p, o.dates, o.values; α=1., lc=:match, lab="observed $(r.label)", color=color_observed, markersize, ylabel, msw=0.)
    if plot_only_fit

        xlo = Dates.value(pd[1]) - 2
        xup = Dates.value(pd[2]) + 2
        # yup = maximum( o.values[pd[1] .< o.dates .< pd[2]] ) * 1.1
        # xlo = Dates.value.(sd) .- 2
        # xup = Dates.value.(ed) .+ 2
        # yup = maximum( o.values[sd .< o.dates .< ed] ) * 1.1
        xlims!(p, xlo, xup)
        # ylims!(p, 0, yup)
    end
    return p
end

# ============================================================================
# latent infections
struct LatentInfectionsPlottingRecipe2{Tsd, Ted, Te, Tl} <: PlottingRecipe
    startdate ::Tsd
    enddate   ::Ted
    expected  ::Te
    label     ::Tl
end

function LatentInfectionsPlottingRecipe2(data::NamedTuple, gp, region, label = "infections")
    e  = expected(data, gp, region, label)
    sd = startdate(data, region)
    ed = enddate(data, region) #ed = data.observations_end


    #isnothing(data.predictors) &&
    return LatentInfectionsPlottingRecipe2( sd, ed, e, label)
    #pstartdate = Date("2020-11-10")
    #return LatentInfectionsPlottingRecipe2( ed, e, p, label, pstartdate)
end

function Plots.plot!(p::Plots.Plot, r::LatentInfectionsPlottingRecipe2; plot_only_fit=false, ylabel = "latent cases", kwargs...)
    e  = r.expected
    ed = r.enddate
    sd = r.startdate

    if !plot_only_fit
        !isnothing(ed) && vline!(p, [ed], lab="end observations", lw=2, lc=:black, hover="$ed")
    end
    plot_confidence_timeseries!(p, e.dates, e.values; label = "latent $(r.label)",  ylabel, kwargs...) #Dict(hover=>strdates)
    if plot_only_fit
        xlo = Dates.value.(sd) .- 1 # xlo = Dates.value(e.dates[1])
        xup = Dates.value.(ed) .+ 1
        # yup = maximum( e.values[e.dates .< ed] ) * 1.1
        xlims!(p, xlo, xup)
        # ylims!(p, 0, yup)
    end
    return p
end

# ============================================================================
# reproduction number for a single region

struct RtPlottingRecipe{Tlo, Tst, Ted, Te, Tl} <: PlottingRecipe
    lockdown  ::Tlo
    startdate ::Tst
    enddate  ::Ted
    expected  ::Te
    label     ::Tl
end

struct RtWithPredictorsPlottingRecipe{Tlo, Tst, Ted, Te, Tl, Psd} <: PlottingRecipe
    lockdown  ::Tlo
    startdate ::Tst
    enddate  ::Ted
    expected  ::Te
    label     ::Tl
    pstartdate::Psd
end

function RtPlottingRecipe(data::NamedTuple, generated_posterior, region, label)
    i  = region isa Integer ? region : findfirst(==(region), data.regions)
    dates      = data.dates[i]
    values     = generated_posterior.rts[i]
    expected   = (; dates, values )
    lo = haskey(data, :lockdown) && !isnothing(data.lockdown) ? Date(data.lockdown) : nothing
    # sd = first(data.dates[i])
    # ed = Date(data.observationsend[i])
    sd = startdate(data, region)
    ed = enddate(data, region)

    isnothing(data.predictors) && return RtPlottingRecipe( lo, sd, ed, expected, label)
    psd = inference_boundaries(data.data[region]) #Date("2020-11-10")
    return RtWithPredictorsPlottingRecipe( lo, sd, ed, expected, label, psd)
end

function Plots.plot!(p::Plots.Plot, r::RtPlottingRecipe; plot_only_fit = false,  kwargs...)
    e  = r.expected
    sd = r.startdate
    ed = r.enddate
    lo = r.lockdown

    if !plot_only_fit
        !isnothing(ed) && vline!(p, [ed], lab="end observations", lw=2, lc=:black, hover="$ed")
        !isnothing(lo) && vline!(p, [lo], lab="lockdown", lw=2, lc=:black, hover="$lo", ls=:dash)
    end
    plot_confidence_timeseries!(p, e.dates, e.values; label = r.label,  kwargs...) #Dict(hover=>strdates)
    if plot_only_fit
        xlo = Dates.value.(sd)
        xup = Dates.value.(ed)
        # yup = maximum( o.values[o.dates .< ed] ) * 1.1
        xlims!(p, xlo, xup)
        # ylims!(p, 0, yup)
    end
    return p
end

function Plots.plot!(p::Plots.Plot, r::RtWithPredictorsPlottingRecipe; plot_only_fit = false,  kwargs...)
    e  = r.expected
    sd = r.startdate
    ed = r.enddate
    lo = r.lockdown
    pd = r.pstartdate

    if !plot_only_fit
        vline!(p, [ed], lab="end observations", lw=2, lc=:black, hover="$ed")
        !isnothing(lo) && vline!(p, [lo], lab="lockdown", lw=2, lc=:black, hover="$lo", ls=:dash)
        vspan!(p, pd, label="fit to survey data", α=0.2, fc=:lightgray)
    end
    plot_confidence_timeseries!(p, e.dates, e.values; label = r.label,  kwargs...) #Dict(hover=>strdates)
    if plot_only_fit
        xlo = Dates.value(pd[1]) - 2
        xup = Dates.value(pd[2]) + 2
        # yup = maximum( o.values[pd[1] .< o.dates .< pd[2]] ) * 1.1
        # xlo = Dates.value.(sd) .- 2
        # xup = Dates.value.(ed) .+ 2
        # yup = maximum( o.values[sd .< o.dates .< ed] ) * 1.1
        xlims!(p, xlo, xup)
        # ylims!(p, 0, yup)
    end
    return p
end
# ============================================================================
# region plot

struct RegionPlottingRecipe{Tr,Tt,T} <: PlottingRecipe
    recipes::Tr
    titles ::Tt
    region ::T
end

posterior2recipe = OrderedDict(
    :expecteds => ObservationPlottingRecipe,
    :rts       => RtPlottingRecipe,
)

posterior2label = OrderedDict(
    :expecteds      => "expected",
    :rts            => "reproduction number",
    :newly_infected => "newly infected"
)

function RegionPlottingRecipe(data::NamedTuple, generated_posterior, region)
    ks = keys(generated_posterior)
    recipes = Vector{PlottingRecipe}()
    titles  = Vector{String}()
    for (k,recipe) in posterior2recipe
        if k in ks
            label = posterior2label[k]
            r     = recipe(data, generated_posterior, region, label)
            title = label * " for the $region"
            push!(recipes, r)
            push!(titles, title)
        end
    end
    RegionPlottingRecipe(recipes, titles, region)
end

function Plots.plot(r::RegionPlottingRecipe; kwargs...)
    plots  = Vector{Plots.Plot}()
    nplots = length(r.recipes)
    for (recipe, title) in zip(r.recipes, r.titles)
        p = plot(; xaxis = true, legend = :outertopright, title)
        plot!(p, recipe; kwargs...)
        push!(plots, p)
    end
    plot(plots..., layout=(nplots,1), size=(1000, nplots*250), sharex=true, link=:x)
end

# Plots.plot(data::NamedTuple, generated_posterior) =
#     plot( RegionPlottingRecipe(data, generated_posterior) )
# Plots.plot(data::NamedTuple, generated_posterior) = plot!(plot(), data, generated_posterior)

# ============================================================================
# reproduction number for all regions

struct RtsPlottingRecipe{Tlo, Tsd, Ted, Tre, Te} <: PlottingRecipe
    lockdown  ::Tlo
    startdates::Tsd
    enddate   ::Ted
    regions   ::Tre
    expecteds ::Te
end

function RtsPlottingRecipe(data::NamedTuple, generated_posterior)
    dates      = data.dates
    samples    = generated_posterior.rts
    expecteds  = [(dates = d, values = g) for (d,g) in zip(dates, samples)]
    lockdown   = haskey(data, :lockdown) && !isnothing(data.lockdown) ? Date(data.lockdown) : nothing # NOTE use data["lockdown"] instead
    startdates = first.(data.dates)
    enddate    = Date.(data.observationsend) # NOTE use data["observations_end"] instead
    RtsPlottingRecipe( lockdown, startdates, enddate, data.regions, expecteds)
end

function Plots.plot(r::RtsPlottingRecipe; kwargs...)
    plots  = Vector{Plots.Plot}()
    nplots = length(r.regions)
    for (expected, region, startdate, enddate) in zip(r.expecteds, r.regions, r.startdates, r.enddate)
        recipe = RtPlottingRecipe(r.lockdown, startdate, enddate, expected, "reproduction number")
        p = plot(; xaxis = true, legend = :outertopright, title="$region")
        plot!(p, recipe; kwargs...)
        push!(plots, p)
    end
    plot(plots..., layout=(nplots,1), size=(1000, nplots*200), sharex=true, link=:x)
end

# ============================================================================
# observations for all regions

struct RegionsOverviewPlottingRecipe{Tr,Tt,T} <: PlottingRecipe
    recipes::Tr
    titles ::Tt
    region ::T
end


function RegionsOverviewPlottingRecipe(data::NamedTuple, generated_posterior, posterior)
    # i  = region isa Integer ? region : findfirst(==(region), data.regions)
    recipes = Vector{PlottingRecipe}()
    titles  = Vector{String}()
    for (i,region) in enumerate(data.regions)
        label  = posterior2label[posterior]
        recipe = posterior2recipe[posterior]
        r      = recipe(data, generated_posterior, region, label)
        title  = i == 1 ? label*" for the $region" : "$region"
        push!(recipes, r)
        push!(titles, title)
    end
    return RegionsOverviewPlottingRecipe(recipes, titles, data.regions)
end

function Plots.plot(r::RegionsOverviewPlottingRecipe; kwargs...)
    # plots  = Vector{Plots.Plot}()
    nplots = length(r.recipes)

    title = r.titles[1]
    p = plot(; xaxis = true, legend = :outertopright, title)
    plot!(p, r.recipes[1]; kwargs...)
    plots = [p]

    for (recipe, title) in Iterators.drop(zip(r.recipes, r.titles),1)
        p = plot(; xaxis = true, legend = nothing, title)
        plot!(p, recipe; kwargs...)
        push!(plots, p)
    end
    plot(plots..., layout=(nplots,1), size=(1000, nplots*250), sharex=true, link=:x)
end
