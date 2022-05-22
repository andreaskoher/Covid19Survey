above(x, v::Int) = sum(skipmissing(x) .> v) / numnotmissing(x) * 100
function above(x, p::Float64)
    @assert 0 <= p <= 1
    v = Int( quantile(x, p) )
    above(x, v)
end

numnotmissing(x) = sum(x->!ismissing(x), x)

# =============================================================================
# remove outliers

function maskoutliers(data, p::Real=99.9)
    outliers = zeros(size(data,1))
    for c in datacols(data)
        threshold = p isa Number ? percentile(data[:,c], p) : p[c]
        @show "threshold for $c: $threshold"
        outliers .+= (threshold .< data[!,c])
        outliers .+= (data[!,c] .< 0)
    end
    return outliers .== 0
end

function maskoutliers(data, threshold::Union{NamedTuple, Dict})
    outliers = zeros(size(data,1))
    for (c, thr) in pairs(threshold)
        @show "threshold for $c: $thr"
        outliers .+= (thr .< data[!,c])
        outliers .+= (data[!,c] .< 0)
    end
    return outliers .== 0
end

function filteroutliers(data; outlier_threshold=99.9, kwargs...)
    mask = maskoutliers(data, outlier_threshold)
    @show "remove $(sum(mask .== 0)/size(data,1)*100) percent"
    data = data[mask,:]
    return data
end

# =============================================================================
# rolling mean with weights

function rollingmean(xs::AbstractVector, window, weights = ones(length(xs)))
    n  = length(xs)-window+1
    ys = Vector{Float64}(undef, n)
    for i in 1:n
        x, w = xs[i:i+window-1], weights[i:i+window-1]
        mask = @. !ismissing(x)
        @assert sum( mask ) > 3 "nr data points: $(sum( mask )) at i=$i"
        ys[i] = mean(x[mask], FrequencyWeights(w[mask]))
        @assert isfinite(ys[i]) "is not finite at i=$i: $(ys[i]) "
    end
    return ys
end

function rollingmean(df, window; exclude=[:date])
    @assert isodd(window)
    Δ = window ÷ 2
    # smoothed = DataFrame( :date => df.date[1+Δ:end-Δ] )
    smoothed = DataFrame( [c => df[1+Δ:end-Δ, c] for c in exclude] )
    for col in names(df, Not(exclude))
        @info col
        smoothed[!,col] = rollingmean(df[:,col], window, df.nrow)#
        @assert all( isfinite.(smoothed[!,col]) ) "$(findall(x-> !isfinite(x), smoothed[!,col]))"
    end
    disallowmissing(smoothed)
end

function total_contacts!(df)
    df[!,:total] .= 0
    for c in [:friends, :strangers, :colleagues, :family]
        df[!,:total] += df[:,c]
    end
    df
end

function rawcontacts(survey::DataFrame, region; kwargs...)
    @_ survey |>
        filter(:region => ==(region), __) |>
        filteroutliers(__; kwargs...) |>
        total_contacts!
end

function rawcontacts(survey::DataFrame; kwargs...)
    @_ survey |>
        filteroutliers(__; kwargs...) |>
        total_contacts!
end

function rawcontacts(
    ; fname = projectdir( "data/survey/selected_survey_responses.csv")
    , kwargs...
)
    @_ fname |>
        CSV.read(__, DataFrame ) |>
        rawcontacts(__; kwargs...)
end
# ===========================================================================
# main

function survey_predictors(
    ; fname = projectdir( "data/survey/selected_survey_responses.csv")
    , kwargs...
)
    survey = CSV.read( fname, DataFrame )
    survey_regional = Dict()
    for region in unique(survey.region)
        @info "load data region $region"
        regional = rawcontacts(survey, region; kwargs...)

        data_cols = names(regional, Not(["date", "region"]))
        contact_cols = [:total, :strangers, :friends, :family, :colleagues]
        survey_regional[region] = @_ regional |>
             groupby(__, :date) |>
             combine(__,
                :region => ( x -> only(unique(x)) ) => :region,
                data_cols .=>( x->mean(skipmissing(x)) ).=> data_cols,
                contact_cols .=>( x->above(x, 0) ).=> ["$(c)-above0" for c in contact_cols],
                contact_cols .=>( x->above(x, 1) ).=> ["$(c)-above1" for c in contact_cols],
                contact_cols .=>( x->above(x, 2) ).=> ["$(c)-above2" for c in contact_cols],
                contact_cols .=>( x->above(x, 3) ).=> ["$(c)-above3" for c in contact_cols],
                contact_cols .=>( x->above(x, 4) ).=> ["$(c)-above4" for c in contact_cols],
                contact_cols .=>( x->above(x, 5) ).=> ["$(c)-above5" for c in contact_cols],
                contact_cols .=>( x->above(x, 6) ).=> ["$(c)-above6" for c in contact_cols],
                contact_cols .=>( x->above(x, 7) ).=> ["$(c)-above7" for c in contact_cols],
                contact_cols .=>( x->above(x, 8) ).=> ["$(c)-above8" for c in contact_cols],
                contact_cols .=>( x->above(x, 9) ).=> ["$(c)-above9" for c in contact_cols],
                contact_cols .=>( x->above(x, 10) ).=> ["$(c)-above10" for c in contact_cols],#12   #15  #19   #25
                contact_cols .=>( x->above(x, 11) ).=> ["$(c)-above11" for c in contact_cols],
                contact_cols .=>( x->above(x, 12) ).=> ["$(c)-above12" for c in contact_cols],
                contact_cols .=>( x->above(x, 13) ).=> ["$(c)-above13" for c in contact_cols],
                contact_cols .=>( x->above(x, 14) ).=> ["$(c)-above14" for c in contact_cols],
                contact_cols .=>( x->above(x, 15) ).=> ["$(c)-above15" for c in contact_cols],
                contact_cols .=>( x->above(x, 16) ).=> ["$(c)-above16" for c in contact_cols],
                contact_cols .=>( x->above(x, 17) ).=> ["$(c)-above17" for c in contact_cols],
                contact_cols .=>( x->above(x, 18) ).=> ["$(c)-above18" for c in contact_cols],
                contact_cols .=>( x->above(x, 19) ).=> ["$(c)-above19" for c in contact_cols],
                contact_cols .=>( x->above(x, 20) ).=> ["$(c)-above20" for c in contact_cols],
                contact_cols .=>( x->above(x, 21) ).=> ["$(c)-above21" for c in contact_cols],
                contact_cols .=>( x->above(x, 22) ).=> ["$(c)-above22" for c in contact_cols],
                contact_cols .=>( x->above(x, 23) ).=> ["$(c)-above23" for c in contact_cols],
                contact_cols .=>( x->above(x, 24) ).=> ["$(c)-above24" for c in contact_cols],
                contact_cols .=>( x->above(x, 25) ).=> ["$(c)-above25" for c in contact_cols],
                nrow) |>
             sort(__, :date) |>
             rollingmean(__, 7, exclude=[:date, :region, :nrow]) |>
             disallowmissing
    end
    @_ survey_regional |>
        vcat( values(__)... ) |>
        sort(__, [:region, :date])
end
