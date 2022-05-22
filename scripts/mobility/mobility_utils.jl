function normalize!(xs::AbstractVector; kwargs...)
    xs[1] == 0. && (@warn "predictor data stream starts with 0, so we assume it is already normalized."; return nothing)
    xs ./= abs(xs[1])
    # xs .-= xs[1]
    xs
end

function replacemissing(xs)
    for i in eachindex(xs)
        if ismissing(xs[i])
            xs[i] = @_ vcat( xs[i-2:i-1], xs[i+1:i+2] ) |>
                skipmissing |>
                mean
        end
    end
    convert(Vector{Float64}, xs)
end

function custom_mean(xs)
    n = 0
    s = 0.
    for x in xs
        if !ismissing(x)
            n += 1
            s += x
        end
    end
    if n == 0
        @warn "all missing"
        return missing
    else
        m = s / n
        @assert isfinite(m)
        return m
    end
end

function kommune2region(
      old2new_names
    ; google_fname = projectdir("data/mobility/google/Global_Mobility_Report.csv")
)
    eng2dan = Dict(
            "Brondby" => "Brøndby",
            "Copenhagen" => "København",
            "Lyngby Taarbæk" => "Lyngby Taarbæk",
            "Nordfyn" => "Nordfyns",
            "Aarhus" => "Århus"
        )

    @_ google_fname |>
        CSV.File |>
        DataFrame |>
        filter(:country_region => ==("Denmark"), __) |>
        filter( :sub_region_1 => x->!ismissing(x), __) |>
        filter( :sub_region_2 => x->!ismissing(x), __) |>
        DataFrames.select(__, :sub_region_1 => :region, :sub_region_2 => :kommune) |>
        disallowmissing |>
        groupby(__, :kommune) |>
        combine(__, unique) |>
        DataFrames.select(__,
            :region => ByRow(x -> old2new_names[x]) => :region,
            :kommune => ByRow( x->string(first(split(x, " Municipality"))) ) => :kommune,
        ) |>
        DataFrames.transform(__, :kommune => ByRow( x-> x ∈ keys(eng2dan) ? eng2dan[x] : x ) => :kommune ) |>
        DataFrames.transform(__, :kommune => ByRow( x->replace(x, "-"=>" ") ) => :kommune) |>
        append!(__,
            DataFrame( Dict(
            "kommune" => ["Ærø", "Samsø", "Fanø", "Læsø"],
            "region"  => ["south", "central", "south", "north"]
        ))) |>
        Dict([row.kommune => row.region for row in eachrow(__)])
end

function json2long(raw; dateformat = dateformat"yyyy-mm-dd HH:MM:SS")
    num_kommunes = length(raw["locations"])
    data = Array(hcat(raw["data"]...)')
    source = string.(repeat(raw["locations"], inner=num_kommunes))
    target = string.(repeat(raw["locations"], outer=num_kommunes))
    dates = Date.(raw["dates"], dateformat)


    df = DataFrame(
        "date" => Vector{Date}(),
        "source" => Vector{String}(),
        "target" => Vector{String}(),
        "flow" =>  Vector{Float64}()
        )

    for (d, flow) in zip(dates, eachcol(data))
        append!(df, DataFrame( (; date = fill(d, num_kommunes^2), source, target, flow) ) )
    end

    return df
end

function rollingmean(xs::AbstractVector, window)
    n  = length(xs)-window+1
    ys = Vector{Float64}(undef, n)
    for i in 1:n
        x = xs[i:i+window-1]
        @assert sum( ismissing.(x) ) < 4 "nr of missing: $(sum( ismissing.(x) )) at i=$i"
        ys[i] = mean(skipmissing(x))
        @assert isfinite(ys[i]) "infinite at i=$i"
    end
    return ys
end

numericcols(df::DataFrame) =
    filter(x->nonmissingtype( eltype( df[!,x] ) ) <: Number, names(df))


function postprocessing(df, window)
    @assert isodd(window)
    Δ = window ÷ 2
    smoothed = DataFrame( :date => df.date[1+Δ:end-Δ] )
    for col in numericcols(df)
        @info col
        smoothed[!,col] = rollingmean(df[:,col], window)
        smoothed[!,col] ./= abs(smoothed[1,col])
        # smoothed[!,col] .-= 1
        # smoothed[!,col] ./= StatsBase.std(smoothed[!,col])
        if !all( isfinite.(smoothed[!,col]) )
            @error "is not finite at $col:"
            println("$(findall(x-> !isfinite(x), smoothed[!,col]))")
        end
    end
    smoothed
end

function assertconsistency(df)
    @assert all( df.date .==  df.date[1]:Day(1):df.date[end] )
    @assert all( [!any( @. isnan(col) | ismissing(col) | ! isfinite(col) ) for col in eachcol(df) if eltype(col) <: Number] )
    return df
end

function drop_empy_columns(df)
    drop = []
    for n in names(df)
        col = df[!,n]
        all(ismissing.(col)) && push!(drop, n)
    end
    return DataFrames.select(df, Not(drop))
end

function postprocessing_mobilty(mobil)
    regionalmobil = Dict()
    for region in unique(mobil.region)
        @info region
        regional = @_ mobil |>
            filter(:region => ==(region), __) |>
            sort(__, :date) |>
            drop_empy_columns |>
            postprocessing(__, 7) |>
            DataFrames.select(__, All(), AsTable([:google,:apple]) => ByRow(mean) => :applegoogle) |>
            disallowmissing |>
            assertconsistency

        regionalmobil[region] = regional
    end
    longformat = vcat( [constant(data, "region", region) for (region, data) in pairs(regionalmobil)]... )
    sort(longformat, [:region, :date])
end

function constant(df::AbstractDataFrame, col, val)
    df[!, col] .= val
    return df
end
