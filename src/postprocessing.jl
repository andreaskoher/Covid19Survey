# ============================================================================
# post-processing
function save_results(fdir, ps, data, chain)
    ignores = [k for (k,v) in pairs(ps) if ((v isa String) && isempty(v))]
    mkpath(fdir)
    @info "Saving at: $fdir"

    let
        dic = Dict( zip( keys(ps), values(ps) ) )
        fname = normpath(fdir, savename("PARAMS", ps, "csv"; ignores) )
        safesave( fname, DataFrame( dic ) )
        bson( normpath(fdir, "params.bson") ,  dic )
    end
    let
        fname = normpath( fdir, savename("DATA", ps, "bson"; ignores) )
        bson( fname, Dict("data"=>data) )
    end

    fname = normpath( fdir, savename("CHAIN", ps, "jls"; ignores) )
    safesave( fname, chain )

    # try #save source code
    #     Base.Filesystem.cp(projectdir("src"), normpath(fdir, "src"))
    #     Base.Filesystem.cp(projectdir("scripts/run_regional.jl"), normpath(fdir, "run_regional.jl"))
    # catch e
    #     println(e)
    # end

    fname
end

function groupregions(gq::NamedTuple)
    ks = keys(gq)
    vs = values(gq)
    n  = length(first(first(gq)))
    d  = Dict()
    for (k,v) in zip(ks,vs)
        d[k] = [[v[j][i] for j in eachindex(v)] for i in 1:n]
    end
    (; d... )
end

function posterior(model, turing_data, chain)
    m_pred = model(turing_data, :prediction)
    gq = Turing.generated_quantities(m_pred, chain)
    gp = vectup2tupvec(reshape(gq, length(gq)))
    groupregions(gp)
end

function rtstats(Rts, dates; startdate=Date("2020-05-15"), stopdate=nothing)
    ib = findfirst(==(startdate), dates)
    ie = isnothing(stopdate) ? length(dates) : findfirst(==(Date(stopdate)), dates)
    Rt_array = hcat(Rts...)[ib:ie,:]
    qs = [quantile(v, [0.025, 0.25, 0.5, 0.75, 0.975]) for v in eachrow(Rt_array)]
    llq, lq, mq, uq, uuq = (eachrow(hcat(qs...))..., )

    return DataFrame((; dates = dates[ib:ie], llq, lq, mq, uq, uuq))
end

# ============================================================================
# post-processing

struct PostProcessing{F,P,I,C,D,M}
    fdir::F
    ps::P
    ignores::I
    chain::C
    data::D
    model::M
end

function PostProcessing(fdir, ps, ignores, exclude, fname)
    chain = let
        chain = read(fname, Chains)
        not(chain, exclude)
    end
    data = let
        fn_data  = replace(fname, "CHAIN_"=>"DATA_")
        fn_data  = replace(fn_data, ".jls"=>".bson")
        BSON.load(fn_data)["data"]
    end

    PostProcessing(fdir, ps, ignores, chain, data, ps.model)
end

function not(c::Chains, exclude::AbstractVector=[])
    isempty(exclude) && return c
    n = size(c,3)
    s = filter(x->x ∉ exclude, 1:n)
    return c[:,:,s]
end

function parse_fname(fname; warmup = nothing)
    ignores = []
    fdir_raw, parsed_args, _ = parse_savename(fname, parsetypes = (Int, Float64, Bool))
    !isnothing(warmup) && (parsed_args["warmup"] = warmup)
    !("predictors" in keys(parsed_args)) && (parsed_args["predictors"] = nothing; push!(ignores, :predictors))
    parsed_args["model"] = getfield(Covid19Survey, Symbol(parsed_args["model"]) )

    ps = NamedTuple(Symbol(k)=>v for (k,v) in pairs(parsed_args))

    folders = split(fdir_raw, "/")
    fdir = normpath("/", folders[1:end-1]...)

    @info "save post-processing in $fdir"
    return (; fdir, ps, ignores)
end

function savechain(p::PostProcessing)
    @unpack fdir, ps, ignores, chain = p
    fname = normpath( fdir, savename("CHAIN", ps, "jls"; ignores) )
    safesave( fname, chain )
end

function plot_chains(p::PostProcessing; plot_results = false)
    @unpack fdir, ps, ignores, chain = p
    n = filter( x->!occursin(r"latent_Rts", x), String.(names(chain)))
    p = plot(chain[n]);
    fname = normpath( fdir, savename("FIG-CHAINSPLOT", ps, "html"; ignores) )
    savefig(p, fname )
    plot_results && display(p) #(display(p); run(`firefox $(fname)`, wait=false))
    return nothing
end

function plot_means(p::PostProcessing; plot_results=false)
    @unpack fdir, ps, ignores, chain = p
    n = filter( x->!occursin(r"latent_Rts", x), String.(names(chain)))
    p = meanplot(chain[n]);
    fname = normpath( fdir, savename("FIG-MEANPLOT", ps, "html"; ignores) )
    savefig(p, fname )
    plot_results && display(p) #(display(p); run(`firefox $(fname)`, wait=false))
    return nothing
end

function skip_warmup(p::PostProcessing)
    @unpack fdir, ps, ignores, chain, data, model = p
    chain = chain[ps.warmup+1:end,:,:]
    PostProcessing(fdir, ps, ignores, chain, data, model)
end

function diagnostics(p::PostProcessing, plot_results = false)
    @unpack fdir, ps, ignores, chain = p

    rhat = describe(chain) |> first |> DataFrame

    highlighters = (
        Highlighter((data, i, j) -> (j == 7) && abs.(data[i, 7]) > 1.1, bold = true, foreground = :red),
    )
    pretty_table(rhat; highlighters, crop=:none)

    highlighters = (
        HTMLHighlighter((data, i, j) -> (j == 7) && abs.(data[i, 7]) > 1.1, HTMLDecoration(color = "red")),
    )
    table = pretty_table(String, rhat; highlighters, backend = Val(:html) )
    fname = normpath( fdir, savename("DIAGNOSTICS", ps, "html"; ignores) )
    write( fname, table)
    plot_results && display(p) #(display(p); run(`firefox $(fname)`, wait=false))
    println("median Rhat: $(median( rhat[:,7] ))")
    println("max Rhat: $(median( rhat[:,7] ))")

    if ps.chains > 1
        @info "gelman diagnostics"
        gelman = gelmandiag(chain) |> DataFrame

        highlighters = (
            Highlighter((data, i, j) -> (j == 2) && data[i, 2] > 1.2, bold = true, foreground = :red),
        )
        pretty_table(gelman; highlighters, crop=:none)

        highlighters = (
            HTMLHighlighter((data, i, j) -> (j == 2) && data[i, 2] > 1.2, HTMLDecoration(color = "red")),
        )
        table = pretty_table(String, gelman; backend = Val(:html), highlighters )
        fname = normpath( fdir, savename("GELMANDIAG", ps, "html"; ignores) )
        write( fname, table)
        plot_results && display(p) #(display(p); run(`firefox $(fname)`, wait=false))
        println("max gelmandiag: $(maximum( gelman[:,2] ))")
        println("median gelmandiag: $(median( gelman[:,2] ))")
        return (; rhat, gelman)
    else
        return (;rhat)
    end
end

function generate_posterior(p::PostProcessing)
    @unpack fdir, ps, ignores, chain, model, data = p
    fname = normpath( fdir, savename("GENERATED-QUANTITIES", ps, "bson"; ignores) )
    # if isfile(fname)
    #     return BSON.load(fname) |> NamedTuple
    # else
    chains_params = Turing.MCMCChains.get_sections(chain, :parameters)
    generated_posterior = Covid19Survey.posterior(model, data.turing_data, chains_params)
    fname = normpath( fdir, savename("GENERATED-QUANTITIES", ps, "bson"; ignores) )
    dic = Dict( zip( keys(generated_posterior), values(generated_posterior) ) )
    bson( fname ,  dic )
    return generated_posterior
    # end
end

function plot_regions(p, gp; plot_results = false, kwargs...)
    @unpack fdir, ps, ignores, data = p
    for r in data.regions
        recipe = Covid19Survey.RegionPlottingRecipe(data, gp, r)
        p = plot(recipe; kwargs...)
        fname = normpath( fdir, savename("FIG-PREDICTION-$(uppercase(r))", ps, "html"; ignores) )
        savefig( p, fname )
        plot_results && display(p) #(display(p); run(`firefox $(fname)`, wait=false))
    end
end

function plot_rt(p, gp; plot_results = false, kwargs...)
    @unpack fdir, ps, ignores, chain, data = p
    recipe = Covid19Survey.RtsPlottingRecipe(data, gp)
    p = plot(recipe; kwargs...)
    fname = normpath( fdir, savename("FIG-RT", ps, "html"; ignores) )
    savefig( p, fname )
    plot_results && display(p) #(display(p); run(`firefox $(fname)`, wait=false))
    return nothing
end

function save_rt(p, gp)
    @unpack fdir, ps, ignores, chain, data = p
    for (i,r) in enumerate(data.regions)
        df = Covid19Survey.rtstats(
            gp.rts[i],
            data.dates[i];
            startdate= first(data.dates[i]),
            stopdate=data.observationsend[i]
        )
        fname = normpath( fdir, savename("Rt-$(uppercase(r))", ps, "csv"; ignores) )
        save(fname, df)
    end
end

function predictive(turing_data)
    turing_data_dct = convert(Dict, turing_data)
    turing_data_dct = convert(Dict{Symbol, Any}, turing_data_dct)
    turing_data_dct[:cases]   = [similar(h, Missing) for h in turing_data_dct[:cases]]
    turing_data_dct[:hospits] = [similar(h, Missing) for h in turing_data_dct[:hospits]]
    turing_data_dct[:deaths]  = [similar(h, Missing) for h in turing_data_dct[:deaths]]
    return namedtuple(turing_data_dct)
end

# function ArviZ.plot_autocorr(p::PostProcessing; plot_results = false)
#     @unpack fdir, ps, ignores, chain = p
#     c = MCMCChains.get_sections(chain, :parameters)
#     n = filter(
#         x-> !occursin("latent_Rts", x) &&
#             !occursin("effect", x)     &&
#             !occursin("ys", x)     &&
#             !occursin("R0", x) &&
#             !occursin("R1", x), String.(names(c))
#     )
#     ArviZ.plot_autocorr(c; var_names=n);
#     fname = normpath( fdir, savename("FIG-AUTOCORR", ps, "png"; ignores) )
#     plt.gcf().savefig(fname)
#     plot_results && display(p) #(display(p); run(`firefox $(fname)`, wait=false))
#     return nothing
# end

# function ArviZ.plot_pair(p::PostProcessing; plot_results = false)
#     @unpack fdir, ps, ignores, chain = p
#     c = MCMCChains.get_sections(chain, :parameters)
#     n = filter(
#         x-> !occursin("latent_Rts", x) &&
#             !occursin("effect", x)     &&
#             !occursin("ys", x)     &&
#             !occursin("R0", x) &&
#             !occursin("R1", x), String.(names(c))
#     )
#     idata = from_mcmcchains(
#         chain;
#         #coords=Dict("school" => schools),
#         #dims=Dict("y" => ["school"], "σ" => ["school"], "θ" => ["school"]),
#         library="Turing",
#     )
#     plot_pair(
#         idata;
#         #coords=Dict("school" => ["Choate", "Deerfield", "Phillips Andover"]),
#         divergences=true,
#     );
#     fname = normpath( fdir, savename("FIG-PAIRS", ps, "png"; ignores) )
#     plt.gcf().savefig(fname)
#     plot_results && display(p) #(display(p); run(`firefox $(fname)`, wait=false))
#     return nothing
# end


function loglikelihood_matrix(ℓ)
    ns, nc = size(ℓ)
    nr = length(ℓ[1].loglikelihoods)
    no = length(ℓ[1].loglikelihoods[1])
    data = zeros(nr*no,ns,nc)
    for c in 1:nc
        for s in 1:ns
            data[:,s,c] = vcat(ℓ[s,c].loglikelihoods...)
        end
    end
    data
end

function loglikelihoods(pp::PostProcessing)
    @unpack fdir, ps, ignores, chain, model, data = pp

    ch = Turing.MCMCChains.get_sections(chain, :parameters)
    m  = model(data.turing_data, :loglikelihood)
    ℓ  = Turing.generated_quantities(m, ch)

    fname = normpath( fdir, savename("LOGLIKELIHOODS", ps, "bson"; ignores) )
    bson(fname,
        Dict(
            "dimensions" => [:data, :step, :chain],
            "loglikelihoods" => loglikelihood_matrix(ℓ)
            )
        )
    return nothing
end
# =============================================================================
sampled_effects(gp, predictor, region) = getindex.(gp.effects[region], predictor)
sampled_effects(gp, predictor) = getindex.(gp.grouped_effects |> first, predictor)

function relative_Rt_change(x, pp::Covid19Survey.PostProcessing, args...; R0 = 1, kwargs...)
    invf = pp.data.turing_data.invlink
    f    = pp.data.turing_data.link
    # @assert R0 == 1 #otherwise check regional_effect_estimation.jl

    @_ sampled_effects(args...) |>
        f.(invf(float(R0)) .+ __ * x) |>
        @. ( __ / R0 - 1 ) * 100
end

function get_effect_sizes(x, pp, args...; kwargs...)
    if x == 0
        return sampled_effects(args...)
    else
        return relative_Rt_change(x, pp, args...; kwargs...)
    end
end

function effect_quantiles(args...; kwargs...)
    v     = get_effect_sizes(args...; kwargs...)
    ll,uu = hpdi(v; alpha=0.05)
    l,u   = hpdi(v; alpha=0.5)
    m     = median(v)
    return ll, l, m, u, uu
end

function effect_quantiles(pp::PostProcessing, gp; effect_on_Rt=0, kwargs...)
    data = []
    for (i,p) in enumerate(pp.data.predictors)
        for (j,r) in enumerate(pp.data.regions)
            q = effect_quantiles(effect_on_Rt, pp, gp, i, j; kwargs...)
            push!(data, [q..., p, r])
        end
        r = "grouped"
        q = effect_quantiles(effect_on_Rt, pp, gp, i; kwargs...)
        push!(data, [q..., p, r])
    end

    DataFrame(
        lower95 = getindex.(data, 1),
        lower50 = getindex.(data, 2),
        median = getindex.(data, 3),
        upper50 = getindex.(data, 4),
        upper95 = getindex.(data, 5),
        predictor = getindex.(data, 6),
        region = getindex.(data, 7),
    )
end

# ==========================================================================
Base.getindex(df::DataFrame, col::Union{AbstractString, Integer, Symbol}) =
    df[!, col] |> only

function plot_effect!(p, effects; addlabel = false, ylabel = :region, scale = 1.)
    label1 = addlabel ? "95% HPDI" : nothing
    label2 = addlabel ? "50% HPDI" : nothing
    label3 = addlabel ? "median" : nothing

    ll, l, m, u, uu, n = @_ effects |>
        getindex.(Ref(__), [:lower95, :lower50, :median, :upper50, :upper95, ylabel])

    plot!(p, [ll, uu],[n,n], lw=2*scale, c=:midnightblue, lab=label1)
    plot!(p, [l,u],[n,n], lw=6*scale, c=:midnightblue, lab=label2)
    scatter!(p, [m], [n], mc=:white, lab=label3, msw=2, lc=:black, ms=5*scale)
    return p
end

function plot_regional_effect!(p, effects::DataFrame, region; kwargs...)
    @_ effects |>
        filter(_.region == region, __) |>
        plot_effect!(p, __; kwargs...)
end

function plot_regional_effects(i, num_predictors, effects::DataFrame, regions = Covid19Survey.regions; xlabel = "effect size")

    p = plot(;
          xlabel = i == num_predictors ? xlabel : ""
        , legend = i == 1 ? :outertopright : nothing
        , title = effects.predictor |> first
        , bottom_margin = i == num_predictors ? 0mm : 3mm
        , top_margin = i == 1 ? 3mm : 3mm
    )

    for (j,region) in enumerate(regions)
        plot_regional_effect!(
            p, effects, region;
            addlabel = (i == 1) && (j == 1)
        )
    end

    plot_regional_effect!(p, effects, "grouped")
    return p
end

function plot_regional_effects(effects::DataFrame; xlabel = "change in Rt [\\%]")
    ps = Vector{Plots.Plot}()
    predictors = unique(effects.predictor)
    num_predictors = length(predictors)
    for (i, predictor) in enumerate(predictors)
        p = @_ effects |>
            filter(_.predictor == predictor, __) |>
            plot_regional_effects(i, num_predictors, __, effects.region; xlabel)
        push!(ps, p)
    end
    plot(ps...,
        layout=(num_predictors,1),
        size=(500, num_predictors*250),
        link=:x, sharex=true
    )
end


function plot_grouped_effects(effects;  xlabel = ""
    , legend = :outertopright
    , title = "grouped effects", kwargs...
)
    ps = Vector{Plots.Plot}()
    predictors = unique(effects.predictor)
    num_predictors = length(predictors)

    p = plot(; xlabel, legend, title, kwargs...)

    for (i, predictor) in enumerate(predictors)
        @_ effects |>
            filter(_.predictor == predictor, __) |>
            filter(_.region == "grouped", __) |>
            Covid19Survey.plot_effect!(p, __; addlabel = i == 1, ylabel = :predictor)
    end
    return p
end

function plot_effects(pp::PostProcessing, gp; plot_results = true, grouped = false, effect_on_Rt = 0., kwargs...)
    @unpack fdir, ps, ignores, chain = pp
    effects = effect_quantiles(pp, gp; effect_on_Rt, kwargs...)

    suffix = "-EFFECTS"
    suffix = grouped ? suffix *= "-GROUPED" : suffix *= "-REGIONAL"
    suffix = effect_on_Rt == 0 ? suffix *= "-ABSOLUTE" :  suffix *= "-RT=$(effect_on_Rt)"
    xlabel = effect_on_Rt == 0 ? "effect size" : "change in Rt [\\%]"
    plot_effects = grouped ? plot_grouped_effects : plot_regional_effects

    fname = normpath( fdir, savename("FIG"*suffix, ps, "html"; ignores) )
    p = plot_effects(effects; xlabel)
    savefig(p, fname )
    plot_results && display(p) #(display(p); run(`firefox $(fname)`, wait=false))

    fname = normpath( fdir, savename("TABLE"*suffix, ps, "csv"; ignores) )
    CSV.write(fname, effects)
    return nothing
end

function plot_effect_crosscorr(pp::PostProcessing, gp; plot_results = true)
    @unpack fdir, ps, ignores, chain, data = pp
    eff  = group(chain, "grouped_effect")
    data = DataFrame( Array(eff), data.predictors )

    begin
        py"""
        import matplotlib.pyplot as plt
        from scipy.stats import pearsonr
        import seaborn as sns
        import matplotlib.cm as cm

        def reg_coef(x,y,label=None,color=None,**kwargs):
            ax = plt.gca()
            r,p = pearsonr(x,y)

            bbox_props = dict(fc=cm.coolwarm((r+1)/2), lw=0, boxstyle="square,pad=1")
            ax.annotate('r = {:.2f}'.format(r), xy=(0.5,0.5), xycoords='axes fraction', ha='center', bbox=bbox_props, size=20)
            ax.set_axis_off()
        """
    end
    #, ci=68, truncate=False, s = 10, color='k', scatter_kws={'alpha':0.3})
    p = begin
        sns = pyimport("seaborn")
        g = sns.PairGrid(Pandas.DataFrame(data))
        g.map_diag(sns.distplot)
        # g.map_lower(py"myregplot") #, ci=68
        g.map_lower(sns.regplot, scatter_kws=Dict("alpha"=>0.1, "s"=> 1))
        g.map_lower(sns.kdeplot, levels=4, color=".2")
        g.map_upper(py"reg_coef")
        plt.gcf()
    end

    suffix = "-EFFECTS-GROUPED-XCORR"
    fname = normpath( fdir, savename("FIG"*suffix, ps, "png"; ignores) )
    plt.savefig(fname)
    plot_results && display(p) #display(p) #(display(p); run(`firefox $(fname)`, wait=false))
    return nothing
end

function runtime(chain)
    return chain.info.stop_time - chain.info.start_time |>
        mean |>
        round |>
        Dates.Second |>
        Dates.CompoundPeriod |>
        Dates.canonicalize
end

function finish(p::PostProcessing, diag)
    @unpack chain = p
    println("number of divervences: $(sum( chain[:,:numerical_error,:]) )")
    println("median tree depth: $( median(chain[:,:tree_depth,:]) )")
    max_tree_depth = maximum(chain[:,:tree_depth,:])
    println("fraction of max tree depth = $(max_tree_depth): $(
         sum( chain[:,:tree_depth,:] .== max_tree_depth ) / size(chain, 1) / size(chain, 3) * 100
    ) \\%")
    println("max Rhat: $(median( diag.rhat[:,7] ))")
    if haskey(diag, :gelman)
        println("max gelmandiag: $(maximum( diag.gelman[:,2] ))")
    end
    @show Covid19Survey.runtime(chain)
end

function postprocessing(fname; plot_results = false, exclude = [], warmup = nothing)
    ## ==========================================================================
    @info "load data"
    fdir, ps, ignores = parse_fname(fname; warmup)
    p = PostProcessing(fdir, ps, ignores, exclude, fname)
    #savechain(p)
    ## ==========================================================================
    @info "plot chain"
    try
        plot_chains(p; plot_results)
    catch e
        @error e; println(e)
    end
    ## ==========================================================================
    @info "meanplot"
    try
        plot_means(p; plot_results)
    catch e
        @error e; println(e)
    end
    ## ==========================================================================
    @info "perform diagnostics"
    # p = skip_warmup(p)
    diag = diagnostics(p)
    ## ==========================================================================
    @info "make predictions"
    gp = generate_posterior(p)
    ## ==========================================================================
    # @info "plot autocorr"
    # plot_autocorr(p; plot_results)
    ## ==========================================================================
    # @info "plot pairs"
    # plot_pair(p; plot_results)
    ## ==========================================================================
    @info "plot regions"
    try
        plot_regions(p, gp; plot_results)
    catch e
        @error e; println(e)
    end
    ## ==========================================================================
    @info "plot rt"
    try
        plot_rt(p, gp; plot_results)
    catch e
        @error e; println(e)
    end
    ## ==========================================================================
    num_covariates = p.data.turing_data.num_covariates
    if num_covariates > 0
        @info "plot predictors"
        # pgfplotsx()
        # default(titlefontsize = 20, legendfontsize = 18, labelfontsize = 18, guidefontsize = 18, tickfontsize = 12, framestyle = :zerolines, yminorgrid = true)
        R0 = 1.5
        try
            Covid19Survey.plot_effects(p, gp; plot_results, grouped = false, effect_on_Rt = 0., R0)
            Covid19Survey.plot_effects(p, gp; plot_results, grouped = false, effect_on_Rt = 1., R0)
            Covid19Survey.plot_effects(p, gp; plot_results, grouped = false, effect_on_Rt = -0.5, R0)

            Covid19Survey.plot_effects(p, gp; plot_results, grouped = true, effect_on_Rt = 0., R0)
            Covid19Survey.plot_effects(p, gp; plot_results, grouped = true, effect_on_Rt = 1., R0)
            Covid19Survey.plot_effects(p, gp; plot_results, grouped = true, effect_on_Rt = -0.5, R0)
        catch e
            @error e; println(e)
        end
        # plotlyjs()
        # default()
        if num_covariates > 1
            try
                Covid19Survey.plot_effect_crosscorr(p, gp; plot_results)
            catch e
                @error e; println(e)
            end
        end
    end
    # -----------------------------------------------------------------------------
    @info "store reproduction number"
    try
        save_rt(p, gp)
    catch e
        @error e; println(e)
    end
    # -----------------------------------------------------------------------------
    @info "calculate log-likelihoods"
    try
        loglikelihoods(p)
    catch e
        @error e; println(e)
    end
    # -----------------------------------------------------------------------------
    finish(p, diag)
    return nothing
end
