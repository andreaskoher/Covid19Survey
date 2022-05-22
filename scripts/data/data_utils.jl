# =========================================================================
#               utils
# ==========================================================================

function continuous_time!(df)
    @assert all( df.date[1]:Day(1):df.date[end] .== df.date )
    df
end

function consistent!(raw; nonnegative = false, positive = false, exclude=:date)
    df = @_ raw |>
        disallowmissing |>
        sort(__, :date) |>
        continuous_time!

    if nonnegative || positive
        dat = @_ df |>
            select(__, Not(exclude)) |>
            Array(__)
        nonnegative && @assert all( dat .>= 0 )
        positive && @assert all( dat .> 0 )
    end
    df
end

# function Base.filter!(conditions::Vector{Pair{Symbol, T}} where T<:Function, df::AbstractDataFrame)
#     for condition in conditions
#         filter!(condition, df)
#     end
#     df
# end
# # Base.filter!(condition::Pair{Symbol, Function}, df::AbstractDataFrame) =
# #     filter!(condition, df)
# Base.filter!(nothing, df::AbstractDataFrame) = df

function constant!(df::AbstractDataFrame, col::Union{Symbol, AbstractString}, val)
    df[:, col] .= val
    df
end

# =========================================================================
#               read data
# ==========================================================================


function read_survey(region; fname, datecol = :date, shift = 0, kwargs...)
    @_ fname |>
        CSV.read(__, DataFrame) |>
        filter(:region => ==(region), __) |>
        select(__, Not(:region)) |>
        rename!(__, datecol => :date) |>
        consistent!(__, nonnegative = true, kwargs...) |>
        transform(__, :date =>( x -> x + Day(shift) )=> :date)
end

function read_mobil(region; fname, kwargs...)
    @_ fname |>
        CSV.read(__, DataFrame) |>
        filter(:region => ==(region), __) |>
        select(__, Not(:region)) |>
        consistent!(__)
end

# read_predictors(region; predictors = All(), predictors_kwargs...)
# read_mobil(region; predictors_kwargs...)



function read_epidata(region; fname, condition = nothing, kwargs...)
    @_ fname |>
        CSV.read(__, DataFrame) |>
        filter(:region => ==(region), __) |>
        filter!(condition, __) |>
        consistent!(__, nonnegative = true, exclude=[:date, :region])
end

function read_predictors(region; predictors = All(), survey_kwargs, mobil_kwargs, kwargs...)
    @_ outerjoin(
            read_survey(region; survey_kwargs...),
            read_mobil(region; mobil_kwargs...);
            on=:date) |>
        select(__, :date, predictors) |>
        sort(__, :date)
end

# =========================================================================
#               Predictor preprocessing
# ==========================================================================


function relativechange!(xs::AbstractVector; kwargs...)
    xs[1] == 0. && (@warn "predictor data stream starts with 0, so we assume it is already normalized."; return nothing)
    xs ./= abs(xs[1])
    xs .-= xs[1]
    return xs
end

function relativechange!(df::AbstractDataFrame; kwargs...)
    for n in names(df, Not(:date))
        relativechange!(df[!,n], kwargs...)
    end
    df
end

function clean_predictors(pred, epid; predictors=names(pred), condition, kwargs...)
    @_ pred |>
        leftjoin(epid, __; on = :date) |>
        sort(__, :date) |>
        filter!(condition, __) |>
        select(__, names(pred)) |>
        disallowmissing(__) |>
        relativechange!(__)
end

# =========================================================================
#               Random Walk
# ==========================================================================

function stepindex(n, j0 = 1; step = 1, kwargs...)
    index = Vector{Int64}(undef, n)
    j = j0
    for i in 1:n
        index[i] = j
        i % step == 0 && ( j += 1 )
    end
    index
end

function earlyphase!(df, startdate)
    j = findlast( <(startdate), df.date)
    !isnothing(j) && (df[1:j,:rw] .= 1)
end

function latephase!(df)
    i = findlast(x->!ismissing(x), df[!, :rw])
    df[i+1:end, :rw] .= df[i, :rw]
end

function last_step(df)
    i = findlast(x->!ismissing(x), df[!, :rw])
    isnothing(i) && return 0
    return df[i, :rw]
end

function walk!(df, sdate::Date, edate::Date; kwargs...)
    i0 = last_step(df) + 1
    is = findfirst(==(sdate), df.date)
    ie = findfirst(==(edate), df.date)
    n  = ie - is + 1
    df[is:ie, :rw] = stepindex(n, i0; kwargs...)
end

function walk!(df, is::Integer, ie::Integer; kwargs...)
    i0 = last_step(df) + 1
    n  = ie - is + 1
    df[is:ie, :rw] = stepindex(n, i0; kwargs...)
end

function wait!(df, sdate::Date, edate::Date)
    i0 = last_step(df)
    is = findfirst(==(sdate), df.date)
    ie = findfirst(==(edate), df.date)
    df[is:ie, :rw] .= i0
end

function wait!(df, is::Integer, ie::Integer)
    i0 = last_step(df)
    df[is:ie, :rw] .= i0
end

function randomwalk!(data, prednames; semiparam = true, startdate = nothing, enddate = nothing, kwargs...)
    rw = DataFrame(:date => data.date, :rw => Vector{Union{Missing, Int64}}(missing, size(data,1)))
    startdate = isnothing(startdate) ? rw.date[1] : startdate
    enddate   = isnothing(enddate) ? data.date[end] : enddate

    earlyphase!(rw, startdate)
    if semiparam && !isnothing(prednames)
        semiparam_rw_model!(rw, data, prednames, startdate, enddate; kwargs...)
    else
        walk!(rw, startdate, enddate; kwargs...)
    end
    latephase!(rw)
    disallowmissing!(rw)
    continuous_time!(rw)
    @assert all( diff(rw.rw) .∈ Ref([0,1]) )
    rw.rw .-= rw.rw[1] .- 1 # @assert rw[1, :rw] == 1

    @_ leftjoin(data, rw; on=:date) |>
        sort(__, :date)# |>
        # disallowmissing!(__)
end

function helperdf(rw, data, prednames, startdate, enddate)
    istart = findfirst(==(startdate), rw.date)
    iend = findfirst(==(enddate), rw.date)

    tmp = DataFrame(
        :date => startdate:Day(1):enddate,
        :index => istart:iend,
        :nonparam => ismissing.(data[istart:iend, prednames[1]])
    )
    tmp
end

function semiparam_rw_model!(rw, args...; kwargs...)
    # rw[!,:nonparam] = ismissing.(data[:, prednames[1]])
    helper = helperdf(rw, args...)
    while size(helper,1) > 0
        helper = rwupdate!(rw, helper; kwargs...)
    end
    rw
end

function rwupdate!(rw, helper; kwargs...)
    istart = helper.index[1]
    @assert ismissing(rw[istart, :rw])
    @assert istart == 1 || !ismissing(rw[istart-1, :rw])

    if helper.nonparam[1]
        i = findfirst(==(false), helper.nonparam)
        n = isnothing(i) ? size(helper,1) : i - 1
        iend = helper[n, :index]
        walk!(rw, istart, iend; kwargs...)
    else
        i = findfirst(==(true), helper.nonparam)
        n = isnothing(i) ? size(helper,1) : i - 1
        iend = helper[n, :index]
        wait!(rw, istart, iend)
    end
    helper[n+1:end, :]
end



function seedingdata(df, seeding_steps; predictors, kwargs...)
    startdate = df.date[1]
    seeds = DataFrame( :date => startdate - Day(seeding_steps) : Day(1) : startdate - Day(1) )
    seeds[:, :region] .= only(unique(df.region))
    seeds[:, :rw] .= 1
    seeds[:, :population] .= only(unique(df.population))
    seeds[:, :inference] .= false
    predictors = string.(predictors)
    for col in names(df, Not([:date, :region, :rw, :population, :inference]))
        seeds[:,col] .= col in predictors ? 0. : missing
    end
    seeds
end

# =========================================================================
#                seeding period
# ==========================================================================

function seedingdata(df, num_impute; predictors, kwargs...)
    startdate = df.date[1]
    seeds = DataFrame( :date => startdate - Day(num_impute) : Day(1) : startdate - Day(1) )
    seeds[:, :region] .= only(unique(df.region))
    seeds[:, :rw] .= 1
    seeds[:, :population] .= only(unique(df.population))
    seeds[:, :inference] .= false
    predictors = isnothing(predictors) ? [] : string.(predictors)
    for col in names(df, Not([:date, :region, :rw, :population, :inference]))
        seeds[:,col] .= col in predictors ? 0. : missing
    end
    seeds
end

function seeding(df; num_impute = nothing, kwargs...)
    isnothing(num_impute) && (return df)
    seeds = seedingdata(df, num_impute; kwargs...)
    df[:, :rw] .+= 1
    vcat(seeds, df)
end

function clip!(df; clip_before_inference = -1, clip_after_inference = -1)
    is = clip_before_inference >= 0 ? findfirst( ==(true) , df.inference) - clip_before_inference : 1
    ie = clip_after_inference >= 0 ? findlast( ==(true) , df.inference) + clip_after_inference : size(df, 1)
    # df.rw -= df.rw[1] - 1
    df[is:ie, :]
end

function coalesce_predictors(data, predictors)
    inference_end = findall(==(-1), diff(data.inference))
    for i in inference_end
        for p in predictors
            data[i+1:end,p] .+= data[i,p]
        end
    end
    data
end

# =========================================================================
#                prediction period
# ==========================================================================

function predictiondata(df; num_predictions, step, predictors, constant_prediction, kwargs...)
    # @show num_predictions, steps, predictors, constant_prediction
    startdate = df.date[end] + Day(1)
    enddate   = startdate + Day(num_predictions - 1)

    prediction = DataFrame( :date => startdate : Day(1) : enddate )
    prediction[:, :region] .= only(unique(df.region))
    prediction[:, :rw] = if constant_prediction
        fill(df.rw[end], num_predictions)
    else
        stepindex(num_predictions, df.rw[end]; step)
    end
    prediction[:, :population] .= only(unique(df.population))
    prediction[:, :inference] .= false
    predictors = isnothing(predictors) ? [] : string.(predictors)
    for col in names(df, Not([:date, :region, :rw, :population, :inference]))
        prediction[:,col] .= col in predictors ? df[end, col] : missing
    end
    prediction
end

function prediction(df; num_predictions, kwargs...)
    num_predictions == 0 && (return df)
    prediction = predictiondata(df; num_predictions, kwargs...)
    vcat(df, prediction)
end


# =========================================================================
#                main function for data merge
# ==========================================================================
function load_data(;
    regions = ["capital", "zealand", "south", "central", "north"],
    kwargs...
)
    predictors = getindex(getindex(kwargs, :predictor_kwargs), :predictors)
    if isnothing(predictors) || (predictors isa Vector && isempty(predictors))
        byregion(regions, merge_data_without_predictors; kwargs...)
    else
        byregion(regions, merge_data; kwargs...)
    end
end

function merge_data(region::Union{Symbol, AbstractString}; epidata_kwargs, predictor_kwargs, rw_kwargs, clip_kwargs = (;), kwargs...)
    pred = read_predictors(region ; predictor_kwargs... )
    predictors = names(pred, Not(:date))
    epid = read_epidata(region; epidata_kwargs...)
    pred = clean_predictors(pred, epid; predictor_kwargs... )
    df = @_ leftjoin(epid, pred; on=:date) |>
                sort(__, :date) |>
                transform!(__, first(predictors) =>( x -> @. !ismissing(x) )=> :inference) |>
                clip!(__; clip_kwargs...)|>
                randomwalk!(__, predictors; rw_kwargs... ) |>
                coalesce.(__, 0.) |>
                coalesce_predictors(__, predictors) |>
                disallowmissing!(__) |>
                seeding(__; predictors, kwargs...) |>
                prediction(__; predictor_kwargs..., rw_kwargs..., kwargs...) |>
                continuous_time!(__)
end

function merge_data_without_predictors(region::Union{Symbol, AbstractString}; epidata_kwargs, predictor_kwargs, rw_kwargs, kwargs...)
    epid = read_epidata(region; epidata_kwargs...)
    df = @_ epid |>
                sort(__, :date) |>
                constant!(__, :inference, false) |>
                randomwalk!(__, nothing; rw_kwargs... ) |>
                disallowmissing!(__) |>
                seeding(__; predictors = nothing, kwargs...) |>
                prediction(__; predictor_kwargs..., rw_kwargs..., kwargs...) |>
                continuous_time!(__)
end

function byregion(regions::AbstractVector, func; kwargs...)
    region = regions[1]
    @show region
    data = func(regions[1]; kwargs...)
    for region in regions[2:end]
        @show region
        data = @_ func(region; kwargs...) |>
            vcat(data, __)
    end
    data
end

# function merge_data(regions::AbstractVector; kwargs...)
#     @show regions[1]
#     data = merge_data(regions[1]; kwargs...)
#     for region in regions[2:end]
#         @show region
#         data = @_ merge_data(region; kwargs...) |>
#             vcat(data, __)
#     end
#     data
# end

# =========================================================================
#
# ==========================================================================






# function clip!(df; clip_before_inference = -1, clip_after_inference = -1)
#     is = clip_before_inference >= 0 ? findfirst( ==(true) , df.inference) - clip_before_inference : 1
#     ie = clip_after_inference >= 0 ? findlast( ==(true) , df.inference) + clip_after_inference : size(df, 1)
#     # df.rw -= df.rw[1] - 1
#     df[is:ie, :]
# end
#
# function constant!(df::AbstractDataFrame, col::Union{Symbol, AbstractString}, val)
#     df[:, col] .= val
#     df
# end
#
# function load_data(regions::AbstractVector; kwargs...)
#     @show regions[1]
#     data = load_data(regions[1]; kwargs...)
#     for region in regions[2:end]
#         @show region
#         data = @_ load_data(region; kwargs...) |>
#             vcat(data, __)
#     end
#     data
# end
#
# function load_data(country::Union{Symbol, AbstractString}, regions::AbstractVector; kwargs...)
#     @show country
#     data = cn = load_data(country; kwargs...)
#     for region in regions
#         @show region
#         data = @_ load_data(region, cn; kwargs...) |>
#             vcat(data, __)
#     end
#     data
# end
#
function find_boundaries(xs::AbstractVector)
    changepoints = diff(xs)
    first(xs) == false ? pushfirst!(changepoints, 0) : pushfirst!(changepoints, 1)
    starts = findall(==(1), changepoints)
    ends = findall(==(-1), changepoints) .- 1
    length(starts) > length(ends) && push!(ends, length(xs))
    return (;starts, ends)
end
#
function inference_boundaries(df::DataFrame)
    bounds = find_boundaries(df.inference)
    startdates = df.date[bounds.starts]
    enddates = df.date[bounds.ends]
    bounds = [[s,e] for (s,e) in zip(startdates, enddates)]
    bounds = vcat(bounds...)
    return bounds
end
#
# function firefox(p; fname="tmp.html", fdir=normpath(homedir(), ".tmp"))
#     # mkpath(dir)
#     filepath = normpath(fdir, fname)
#     savefig(p, filepath)
#     run(`firefox $(filepath)`, wait=false)
#     return p
# end
