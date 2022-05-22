# ===============================================================================
#                           main function
# ===============================================================================

function turingformat(data;
    predictors        = [],
    observable        = "hospit",
    link              = KLogistic(4.),
    invlink           = KLogit(4.),
    kwargs...
)
    info = Dict{String, Any}()
    info["link"]        = link
    info["invlink"]     = invlink
    info["num_impute"]  = getindex(kwargs, :num_impute)
    info["rwstep"]      = getindex( getindex(kwargs, :rw_kwargs), :step)
    info["predictors"]  = (isnothing(predictors) || isempty(predictors) ) ?
        nothing : string.(predictors)
    info["observable"]  = observable
    regions!(info, data)
    dates!(info; kwargs...)
    info["observed"] = turing_observables(info, observable; kwargs...)
    initially_infected!(info; cases_column = "cases", kwargs...)
    covariates!(info)
    randomwalks!(info, data)
    serial_interval!(info; nbins = 15)
    info["inf2hosp"] = inf2hosp(; nbins = 40)

    (;
        turing_data                 = turingdata(info),
        num_impute                  = info["num_impute"],
        data                        = info["data"],
        observationsend             = info["observationsend"],
        observationsstart           = info["startdates"],
        num_time_steps              = info["num_time_steps"] ,
        num_observations            = info["num_observations"],
        regions                     = info["regions"],
        dates                       = [vec(info["data"][region].date) for region in info["regions"]],
        populations                 = info["populations"],
        predictors                  = info["predictors"],
        observables                 = info["observable"],
        num_predictions             = 0,
    )
end



# =========================================================================
#
# ==========================================================================

function regions!(data::Dict, df::DataFrame, observationsend = nothing)
    # @unpack predictors = data
    data["regions"] = unique(df.region)
    data["num_regions"] = length(data["regions"] )
    data["populations"] = unique(df.population)
    regional = Dict{String, DataFrame}()
    for region in data["regions"]
        r = regional[region] = @_ df |>
            filter(:region => ==(region), __) |>
            select(__, Not(:region))
        @assert all( r.date[1]:Day(1):r.date[end] .== r.date )
        @assert r.rw[1] == 1 && issorted(r.rw)
        # !isnothing(predictors) && disallowmissing!( r[!,predictors] )
    end
    data["data"] = regional
    return nothing
end

function dates!(data; num_predictions = 0, observationsend = nothing, kwargs...)
    @unpack num_impute = data
    startdates = Vector{Date}()
    dates      = Vector{Vector{Date}}()
    for region in data["regions"]
        regional = data["data"][region]
        push!(startdates, regional[num_impute + 1, :date])
        push!(dates, regional.date)
    end
    isnothing(observationsend) && ( observationsend = last.(dates) )
    observationsend isa Date && ( observationsend = fill(observationsend, data["num_regions"]) )

    data["observationsend"] = Date.(observationsend)
    data["startdates"] = startdates
    data["num_with_prediction"] = data["num_time_steps"] = length.(dates)
    data["num_without_prediction"] = data["num_with_prediction"] .- num_predictions
    data["dates"] = dates
    data["dates_turing"] = turing_observables(data, :date; kwargs...)
    data["num_observations"]  = length.(data["dates_turing"])
    return nothing
end

function limit(dates::Vector, s::Union{Int, Date}, e::Union{Int, Date})
    is = s isa Int ? s : findfirst(==(s), dates)
    ie = e isa Int ? e : findfirst(==(e), dates)
    dates[is:ie]
end

function limit(df::DataFrame, s::Union{Int, Date}, e::Union{Int, Date})
    is = s isa Int ? s : findfirst(==(s), df.date)
    ie = e isa Int ? e : findfirst(==(e), df.date)
    df[is:ie, :]
end

function turing_observables(region, data::Dict, observable; skipseed = true, tovec=true, kwargs...)
    @unpack startdates, observationsend, regions = data
    i = findfirst(==(region), regions)
    s = skipseed ? startdates[i] : 1
    e = observationsend[i]
    d = @_ data["data"][region] |>
        limit(__, s, e) |>
        select(__, observable) |>
        disallowmissing |>
        Array
    return tovec ? tryvec(d) : d
end

function turing_observables(region, data::DataFrame, observable; regions, startdates = nothing, observationsend = nothing, skipseed = true, tovec=true, kwargs...)
    regional = filter(:region => ==(region), data)
    i = findfirst(==(region), regions)
    s = skipseed && !isnothing(startdates) ? startdates[i] : 1
    e = !isnothing(observationsend) ? observationsend[i] : size(regional, 1)
    d = @_  regional |>
        limit(__, s, e) |>
        select(__, observable) |>
        disallowmissing |>
        Array
    return tovec ? tryvec(d) : d
end

function tryvec(xs)
    size(xs, 2) == 1 && return vec(xs)
    xs
end

function turing_observables(data, observable; kwargs...)
    # @unpack num_regions = data
    regions = data isa Dict ? data["regions"] : getindex( kwargs, :regions )
    o = turing_observables(regions[1], data, observable; regions, kwargs...)
    obs = [o]
    for region in regions[2:end]
        o = turing_observables(region, data, observable; regions, kwargs...)
        push!(obs, o)
    end
    obs
end


function initially_infected!(info::Dict; cases_column = "cases", kwargs...)
    cases = turing_observables(info, cases_column; kwargs...)
    init_infected = first.( cases )
    init_infected[init_infected .< 10] .= 10
    info["init_infected"] = init_infected
    return nothing
end

function covariates!(info)
    @unpack num_regions, predictors = info
    if isnothing(predictors)
        info["covariates"]       = [Array{Float64,2}(undef, 1,1) for i in 1:num_regions]
        info["num_covariates"]   = 0
    else
        info["num_covariates"]   = length(predictors)
        info["covariates"] = @_ info |>
            turing_observables(__, predictors, skipseed = false, tovec=false) |>
            convert(Vector{Matrix{Float64}}, __)
    end
    return nothing
end

function randomwalks!(data, merged_data)
    rw = Vector{Vector{Int64}}()
    for region in data["regions"]
        push!(rw,
            @_ merged_data |>
                filter(:region => ==(region), __) |>
                __.rw
        )
    end
    data["rt_step_indices"] = rw
    data["num_rt_steps"] = last.(rw)
    data["rwscale"] = sqrt(data["rwstep"])
    return nothing
end

holiday(
    date::Date;
    specialdays = [Date("2020-12-24"), Date("2020-12-25"), Date("2020-12-31"), Date("2021-01-01")]
) = date ∈ specialdays

function serial_interval!(data; nbins = 15)
    data["serial_interval"] = serial_interval(; nbins)
    data["num_si"] = nbins
    return nothing
end;

function serial_interval(; nbins = 15)
    # Imperial Report
    r = vcat([0], 1.5:1:(nbins + .5))
    # diff(cdf.(GammaMeanStd(6.5, .62), r))
    # return p / sum(p)
    ## Brauer et al. Science (2020)
    p = diff(cdf.(GammaMeanStd(5.06, 2.11), r))
    return p# / sum(p)
end

function inf2hosp(; nbins = 40)
    x = rand(GammaMeanCv(5.1, 0.86), 100_000) + rand( Weibull(0.845, 5.506), 100_000)
    hist = fit(Histogram, x,  [0, 1.5:1:(nbins + .5)...], closed=:left)
    hist.weights / sum(hist.weights)
end

function inf2case(nbins = 30)
    x = rand(GammaMeanCv(5.1, 0.86), 100_000)
    hist = fit(Histogram, x,  [0, 1.5:1:(nbins + .5)...], closed=:left)
    hist.weights / sum(hist.weights)
end

# ============================================================================
function turingdata(info)
    vars = ["num_time_steps", "num_observations", "num_regions", "num_rt_steps",
           "rwscale", "rt_step_indices", "num_impute", "observed", "inf2hosp",
           "link", "invlink", "serial_interval", "populations", "dates_turing",
           "covariates", "num_si", "init_infected", "inf2hosp", "num_covariates"]

    d = Dict()
    for v in vars
        d[Symbol(v)] = info[v]
    end
    return (; d...)
end
