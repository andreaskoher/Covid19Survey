function process_ssi(df, value_name; dan2eng, regions, kwargs...)
    @_ df |>
        DataFrames.select(__, 1=>:region, 2=>:date, 3=>:counts) |>
        filter(:region => x->!ismissing(x), __) |>
        transform(__, :region => ( x->replace(x, dan2eng...) )=> :region) |>
        filter(:region => x -> x ∈ regions, __) |>
        unstack(__, :date, :region, :counts) |>
        insert_missing_dates!(__) |>
        continuous_time!(__) |>
        coalesce.(__, 0) |>
        disallowmissing!(__) |>
        DataFrames.stack(__, regions; variable_name = :region, value_name)
end

function insert_missing_dates!(df)
    ts = df.date[1] : Day(1) : df.date[end]
    missing_date = .! in.(ts, Ref(df.date))
    if any(missing_date)
        for t in ts[missing_date]
            @warn "missing date: $t"
            insertrow!(df, t)
        end
    end
    df
end

function continuous_time!(df)
    @assert all( df.date[1]:Day(1):df.date[end] .== df.date )
    df
end

function insertrow!(df, d::Date, v=0)
    i = findfirst(==(d-Day(1)), df.date)
    v = [eltype(c) == Date ? d : v for c in eachcol(df)]
    insert!.(eachcol(df), i+1, v)
    nothing
end

function Base.filter!(conditions::Vector{Pair{Symbol, Function}}, df::AbstractDataFrame)
    for condition in conditions
        filter!(condition, df)
    end
    df
end
Base.filter!(nothing, df::AbstractDataFrame) = df


function replace_initial_missing!(df, v=0)
    for n in ["cases", "hospit", "deaths"]
        i = findfirst(x->!ismissing(x), df[!,n]) - 1
        i > 0 && (df[1:i,n] .= v)
    end
    df
end

function limit!(df, conditions)
    if conditions isa Pair
        filter!(conditions, df)
    elseif conditions isa Vector && !isempty(conditions)
        filter!.(conditions, Ref(df))
    end
    return df
end

function join_epidata(cases, hospit, deaths; kwargs...)
    @_ outerjoin(cases, hospit; on=[:region, :date]) |>
        outerjoin(__, deaths; on=[:region, :date]) |>
        sort(__, :date) |>
        byregion(__, consistent; kwargs...)
end

function epidemicstart!(df, epidemicstart)
    startdate = df.date[1]
    initphase = DataFrame( :date => startdate - Day(epidemicstart) : Day(1) : startdate - Day(1) )
    for col in names(df, Not(:date))
        initphase[!,col] .= missing
    end
    vcat(initphase, df)
end

function consistent(df, region; conditions, kwargs...)
    @_ df |>
        filter(:region => ==(region), __) |>
        replace_initial_missing!(__, 0) |>
        limit!(__, conditions) |>
        continuous_time!(__) |>#
        disallowmissing!(__) |>
        constant!(__, :region, region)
end

function byregion(df, func; kwargs...)
    regions = unique(df.region)
    region = pop!(regions)
    newdf = func(df, region; kwargs... )
    for region in regions
        transformed = @_ func(df, region; kwargs... )
        newdf = vcat(newdf, transformed)
    end
    newdf
end

function smoothing!(df; kwargs...)
    for col in ["cases", "hospit", "deaths"]
        df[:,"$(col)_smoothed"] = rollingmean(df[:,col]; kwargs...)
    end
    df
end

function rollingmean(xs; window = 7, kwargs...)
    @assert isodd(window)
    Δ  = window ÷ 2
    n  = length(xs)
    ys = Vector{Float64}(undef, n)
    for i in 1:n
        is = max(1, i-Δ)
        ie = min(n, i+Δ)
        ys[i] = mean(xs[is:ie])
    end
    return ys
end

function constant!(df::AbstractDataFrame, col::Union{Symbol, String}, val)
    df[:, col] .= val
    df
end

function postprocess(df, region; populations, kwargs...)
    @_ df |>
        filter(:region => ==(region), __) |>
        constant!(__, :population, populations[region]) |>
        postprocess
end

function incidence(cases, pop)
    n = only(unique(pop))
    T = length(cases)
    incidence = zeros(T)
    for t in 1:T
        incidence[t] = mean(cases[ max(1, t-6) : t ]) * 7 / n * 100_000
    end
    incidence
end

function postprocess(df; kwargs...)
    @_ df |>
        smoothing!(__; kwargs...) |>
        transform(__,
            [:cases, :population] => ( (x,n) -> @. x / n * 100_000) => :cases_per_100000,
            [:hospit, :population] => ( (x,n) -> @. x / n * 100_000) => :hospit_per_100000,
            [:deaths, :population] => ( (x,n) -> @. x / n * 100_000) => :deaths_per_100000,
            [:cases, :population] => ( (x,n) -> incidence(x,n)) => :cases_incidence,
            [:hospit, :population] => ( (x,n) -> incidence(x,n)) => :hospit_incidence,
            [:deaths, :population] => ( (x,n) -> incidence(x,n)) => :deaths_incidence,
            :cases  => cumsum => :cumulative_cases,
            :hospit => cumsum => :cumulative_hospit,
            :deaths => cumsum => :cumulative_deaths,
            )
end

function epidata(;
    fname_cases  = projectdir( "data/epidata/08_bekraeftede_tilfaelde_pr_dag_pr_regions.csv" ),
    fname_hospit = projectdir( "data/epidata/06_nye_indlaeggelser_pr_region_pr_dag.csv" ),
    fname_deaths = projectdir( "data/epidata/07_antal_doede_pr_dag_pr_region.csv" ),
    kwargs...
)

    cases = @_ fname_cases |>
        CSV.read(__, DataFrame ) |>
        process_ssi(__, :cases; kwargs...)


    hospit = @_ fname_hospit |>
        CSV.read(__, DataFrame ) |>
        process_ssi(__, :hospit; kwargs...)

    deaths = @_ fname_deaths |>
        CSV.read(__, DataFrame ) |>
        process_ssi(__, :deaths; kwargs...)

    regional = @_ join_epidata(cases, hospit, deaths; kwargs...) |>
        byregion(__, postprocess; kwargs...)

    @_ regional |>
        groupby(__, :date) |>
        combine(__, :region => x->"denmark", :cases => sum, :hospit => sum, :deaths => sum, :population => sum; renamecols=false ) |>
        postprocess |>
        vcat(regional, __) |>
        sort(__, [:region, :date])
end
