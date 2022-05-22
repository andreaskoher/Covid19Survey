@model function parametricmodel(
    θ,
    cntxt::Symbol = :inference, # choose from :loglikelihood, :inference, :prediction
    ::Type{TV} = Vector{Float64},
    ::Type{V}  = Float64
) where {TV, V}

    @unpack num_time_steps,
            num_observations,
            num_regions,
            num_rt_steps,
            invlink,
            num_covariates,
            rwscale,
            init_infected,
            observed = θ

    ############# 1.) time varying reproduction number
    R0s        ~ filldist(truncated(Normal(1., 0.1), 0., 4.), num_regions) # initialize rt around rt = 1
    σ_rt       ~ truncated(Normal(0.3*rwscale, .02*rwscale), 0, .5*rwscale) # see Unwin et al.

    # draw latent random walk
    latent_rts_z = [TV(undef, num_rt_steps[m]) for m in 1:num_regions]
    latent_rts = [TV(undef, num_rt_steps[m]) for m in 1:num_regions]
    for m in 1:num_regions
        μ_rt = invlink(R0s[m])
        latent_rts_z[m] ~ RandomWalk(num_rt_steps[m], σ_rt, μ_rt, invlink)
        latent_rts[m] = latent_rts_z[m] * σ_rt .+ μ_rt
    end

    # draw effect sizes
    grouped_effect ~ filldist( Laplace(0, .2), num_covariates) # national effect sizes
    effect_std ~ GammaMeanStd(0.03, 0.02) # regional variation
    effects_z ~ filldist( MvNormal( num_covariates, 1.), num_regions)
    effects = [ effects_z[:,m] .* effect_std + grouped_effect for m in 1:num_regions] # regional effects

    # construct rt from latent random walk and predictors
    rts = TV[TV(undef, num_time_steps[m]) for m in 1:num_regions]
    reproductionnumber!(rts, latent_rts, effects, θ)
    # CovidSurvey.random_walk_model!(rts, θ, latent_rts)

    ############ 2.) infection dynamics
    ys ~ arraydist(Exponential.(3 * init_infected)) # initial infected
    newly_infecteds = TV[TV(undef, num_time_steps[m]) for m in 1:num_regions]
    infections!(newly_infecteds, ys, rts, θ)

    ########### 4.) expected hospitalizations
    ihr   ~ truncated(Normal(2.8/100, 0.2/100), 0, 5/100)
    ϕ     ~ GammaMeanStd(50, 20)

    expecteds = TV[TV(undef, num_observations[m]) for m in 1:num_regions]
    expected!(expecteds, newly_infecteds, θ, ihr)

    ########### 4.) compare expected to observed hospitalizations
    observation_model = NegativeBinomialModel(θ, ϕ, expecteds)
    if cntxt == :inference
        likelihood = logpdf(observation_model, observed)
        Turing.@addlogprob! likelihood
    elseif cntxt == :loglikelihood
        loglikelihoods = pointwise_loglikelihoods(observation_model, observed)
        (; loglikelihoods)
    elseif cntxt == :prediction
        (
            ; newly_infecteds
            , expecteds
            , rts
            , predicted = prediction(observation_model)
            , grouped_effects = [grouped_effect for m in 1:num_regions]
            , effects
        )
    else
        @error "choose cntxt from [:inference, :loglikelihood, :prediction]"
    end
end

@model function nonparametricmodel(
    θ,
    cntxt::Symbol = :inference, # choose from :loglikelihood, :inference, :prediction
    ::Type{TV} = Vector{Float64},
    ::Type{V}  = Float64
) where {TV, V}

    @unpack num_time_steps,
            num_observations,
            num_regions,
            num_rt_steps,
            invlink,
            rwscale,
            init_infected,
            observed = θ

    ############# 1.) time varying reproduction number
    R0s        ~ filldist(truncated(Normal(1., 0.1), 0., 4.), num_regions) # initialize rt around rt = 1
    σ_rt       ~ truncated(Normal(0.3*rwscale, .02*rwscale), 0, .5*rwscale) # see Unwin et al.

    # draw latent random walk
    latent_rts_z = [TV(undef, num_rt_steps[m]) for m in 1:num_regions]
    latent_rts = [TV(undef, num_rt_steps[m]) for m in 1:num_regions]
    for m in 1:num_regions
        μ_rt = invlink(R0s[m])
        latent_rts_z[m] ~ RandomWalk(num_rt_steps[m], σ_rt, μ_rt, invlink)
        latent_rts[m] = latent_rts_z[m] * σ_rt .+ μ_rt
    end

    # construct rt from latent random walk and predictors
    rts = TV[TV(undef, num_time_steps[m]) for m in 1:num_regions]
    reproductionnumber!(rts, latent_rts, θ)
    # CovidSurvey.random_walk_model!(rts, θ, latent_rts)

    ############ 2.) infection dynamics
    ys ~ arraydist(Exponential.(3 * init_infected)) # initial infected
    newly_infecteds = TV[TV(undef, num_time_steps[m]) for m in 1:num_regions]
    infections!(newly_infecteds, ys, rts, θ)

    ########### 4.) expected hospitalizations
    ihr   ~ truncated(Normal(2.8/100, 0.2/100), 0, 5/100)
    ϕ     ~ GammaMeanStd(50, 20)

    expecteds = TV[TV(undef, num_observations[m]) for m in 1:num_regions]
    expected!(expecteds, newly_infecteds, θ, ihr)

    ########### 4.) compare expected to observed hospitalizations
    observation_model = NegativeBinomialModel(θ, ϕ, expecteds)
    if cntxt == :inference
        likelihood = logpdf(observation_model, observed)
        Turing.@addlogprob! likelihood
    elseif cntxt == :loglikelihood
        loglikelihoods = pointwise_loglikelihoods(observation_model, observed)
        (; loglikelihoods)
    elseif cntxt == :prediction
        (
            ; newly_infecteds
            , expecteds
            , rts
            , predicted = prediction(observation_model)
        )
    else
        @error "choose cntxt from [:inference, :loglikelihood, :prediction]"
    end
end


# ============================================================================
#                   UTILS
# ============================================================================
# random walk model

function reproductionnumber!(rts, latent_rts, θ::NamedTuple)
	@unpack num_regions, rt_step_indices, link = θ
	for m in 1:num_regions
		reproductionnumber!(
			rts[m],
			latent_rts[m],
			rt_step_indices[m],
			link
		)
	end
end

function reproductionnumber!(rt::AbstractVector{<:Real}, latent_rt, rt_step_index, link)
	for i in eachindex(rt)
		rt[i] = link( latent_rt[ rt_step_index[i] ] ) #NOTE @inbounds
	end
end

# ============================================================================
# random-walk model with predictors

function reproductionnumber!(rts, latent_rts, effect_sizes, θ::NamedTuple)
	@unpack num_regions, link, rt_step_indices, covariates = θ #rt_step_indices
	for m in 1:num_regions
		reproductionnumber!(
			rts[m],
			latent_rts[m],
			rt_step_indices[m],
			link,
			covariates[m],
			effect_sizes[m]
		)
	end
	return nothing
end

function reproductionnumber!(rt::AbstractVector{<:Real}, latent_rt, rt_step_index, link, covariate, effect_size)
	for i in eachindex(rt)
		x = latent_rt[ rt_step_index[i] ] #NOTE @inbounds
		for j in eachindex(effect_size)
			x += covariate[i,j] * effect_size[j]
		end #NOTE @inbounds
		rt[i] = link( x ) #NOTE @inbounds
	end
end

# ============================================================================
# susceptible-infected model
function infections!(newly_infecteds, ys, rts, θ)
	@unpack num_regions = θ
	for m in 1:num_regions
		initepidemic!(newly_infecteds[m], rts[m], ys[m], θ)
		runepidemic!(newly_infecteds[m], rts[m], θ)
	end
end

function initepidemic!(newly_infected, rt, y, θ)
	@unpack num_impute = θ
	newly_infected[1] = y
	for t in 2:num_impute
		newly_infected[t] = y
	end
end

function runepidemic!(newly_infected, rt, θ)
	@unpack num_impute, serial_interval, num_si = θ
	num_time_step = length(newly_infected)

	for t = (num_impute + 1):num_time_step
		effectively_infectious = sum(newly_infected[τ] * serial_interval[t - τ] for τ = (t - 1):-1:max(t-num_si,1))
		newly_infected[t] = effectively_infectious * rt[t]
	end
	return nothing
end

# ============================================================================
# observations

function expected!(expecteds::Vector{Vector{T}}, newly_infecteds::Vector{Vector{T}}, θ, args...) where T<:Real
	@unpack num_regions = θ
	for m in 1:num_regions
		expected!(
			expecteds[m],
			newly_infecteds[m],
			θ, args...
		)
	end
	return nothing
end

function expected!(expecteds::AbstractVector{T}, newly_infecteds::AbstractVector{T}, θ, α) where T<:Real
	@unpack num_regions, num_impute, inf2hosp = θ
	delay_length = length(inf2hosp)
	for i in eachindex(expecteds)
	    t = i + num_impute
		expecteds[i] = α * sum(newly_infecteds[τ] * inf2hosp[t - τ] for τ = (t - 1):-1:max(t-delay_length,1))
	end
	return nothing
end

# ============================================================================
# observation model

struct NegativeBinomialModel{T,P,E}
	θ::T
	ϕ::P
	expecteds::E
end

function _logpdf(ys, μs, ϕ)
	dist = arraydist(NegativeBinomial2.(μs, Ref(ϕ)))
	return logpdf(dist, ys)
end

function Distributions.logpdf(obsmodel::NegativeBinomialModel{T,P,E}, observed) where {T,P,E}
	@unpack θ, ϕ, expecteds = obsmodel
	@unpack num_regions = θ

	ℓ = zero(P)
	for m in 1:num_regions
		μs = expecteds[m]
		ys = observed[m]
		ℓ += _logpdf(ys, μs, ϕ)
	end
	return ℓ
end

# ============================================================================
# pointwise_loglikelihoods
function pointwise_loglikelihoods(observation_model::NegativeBinomialModel, observed)
	@unpack θ, ϕ, expecteds = observation_model
	@unpack num_regions, num_observations = θ

	logliks = [Vector{Float64}(undef, num_observations[m]) for m in 1:num_regions]
	for m in 1:num_regions
		for (i, (e,o)) in enumerate( zip(expecteds[m], observed[m]) )
			logliks[m][i] = logpdf(NegativeBinomial2(e, ϕ), o)
		end
	end
	return logliks
end

# ============================================================================
# predictions
function prediction(obsmodel::NegativeBinomialModel)
	@unpack θ, ϕ, expecteds = obsmodel
	@unpack num_regions, num_observations = θ

	predictions = [Vector{Float64}(undef, num_observations[m]) for m in 1:num_regions]
	for m in 1:num_regions
		es = expecteds[m]
		es[es .< 1e-2] .= 1e-2
        es[es .> typemax(Int64)] .= 1e12
		for (i, e) in enumerate(es)
			predictions[m][i] = rand(NegativeBinomial2(e, ϕ))
		end
	end
	return predictions
end
