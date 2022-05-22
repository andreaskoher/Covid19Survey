"""
This is a generalization of the logit-link, which is referred to as the scaled-logit.
That is to say we expect the quantity of interest x (the reproduction number R,
for instance) to have some carry capacity.

```math
\\sigma_K = KLogistic(K)
```

```math
\\sigma_K(x) = \\frac{K}{e^{-x} + 1}
```
"""
struct KLogistic{T}
	k::T
end
(f::KLogistic)(x::T)  where {T} = inv(exp(-x) + one(x)) * T(f.k)

struct KLogit{T}
	k::T
end
(f::KLogit)(x::T)  where {T} = log(x / (T(f.k) - x))

"""
    GammaMeanCv(mean, cv)

Gamma parametrized by the mean μ and coefficient of variation cv = σ/μ.

In terms of Turings parametrisation we have the following relationship:
μ = αθ
cv = 1 / √α

## References
- https://www.rdocumentation.org/packages/EnvStats/versions/2.3.1/topics/GammaAlt
"""
function GammaMeanCv(mean, cv)
    α = cv^(-2)
    θ = mean / α
    return Gamma(α, θ)
end

"""
    GammaMeanStd(mean, std)

Gamma parametrized by the mean μ and standard deviation σ.

In terms of Turings parametrisation we have the following relationship:
μ = αθ
σ² = αθ²
"""
function GammaMeanStd(μ, σ)
    α = (μ/σ)^2
    θ = σ^2/μ
    return Gamma(α, θ)
end

struct RandomWalk{Tn, Ts, Tx, F} <: ContinuousMultivariateDistribution
    n::Tn
    s::Ts
    x0::Tx
	f::F
end
RandomWalk(n) = RandomWalk(n, 1., 0.)

Distributions.rand(rng::AbstractRNG, d::RandomWalk) = begin
    x = Vector{typeof(d.x0)}(undef, d.n)
    Distributions._rand!(rng, d, x)
    return x
end

Distributions._rand!(rng::AbstractRNG, d::RandomWalk, x::AbstractVector) = begin
	x .= ( randn(d.n)*1e-2 .+ d.f(1.) .- d.x0 ) / d.s
	return nothing
end

Distributions.logpdf(
    d::RandomWalk, x::AbstractVector
) = begin
	ℓ  = logpdf( Normal(d.x0, d.s), first(x) )
	ℓ += logpdf( MvNormal( d.n-1, d.s ), diff(x) )
	return ℓ
end
Distributions._logpdf( d::RandomWalk, x::AbstractVector) = logpdf( d, x)

Turing.Bijectors.bijector(d::RandomWalk) = Bijectors.Identity{1}()

Base.length(d::RandomWalk) = d.n


"""
    NegativeBinomial2(μ, ϕ)

Mean-variance parameterization of `NegativeBinomial`.

## Derivation
`NegativeBinomial` from `Distributions.jl` is parameterized following [1]. With the parameterization in [2], we can solve
for `r` (`n` in [1]) and `p` by matching the mean and the variance given in `μ` and `ϕ`.

We have the following two equations

(1) μ = r (1 - p) / p
(2) μ + μ^2 / ϕ = r (1 - p) / p^2

Substituting (1) into the RHS of (2):
  μ + (μ^2 / ϕ) = μ / p
⟹ 1 + (μ / ϕ) = 1 / p
⟹ p = 1 / (1 + μ / ϕ)
⟹ p = (1 / (1 + μ / ϕ)

Then in (1) we have
  μ = r (1 - (1 / 1 + μ / ϕ)) * (1 + μ / ϕ)
⟹ μ = r ((1 + μ / ϕ) - 1)
⟹ r = ϕ

Hence, the resulting map is `(μ, ϕ) ↦ NegativeBinomial(ϕ, 1 / (1 + μ / ϕ))`.

## References
[1] https://reference.wolfram.com/language/ref/NegativeBinomialDistribution.html
[2] https://mc-stan.org/docs/2_20/functions-reference/nbalt.html
"""
function NegativeBinomial2(μ, ϕ)
    p = 1 / (1 + μ / ϕ)
    r = ϕ
	# !(0 < p <= 1) && (println("μ:$(μ.value) ϕ:$(ϕ.value)"))
    return NegativeBinomial(r, p)
end

DrWatson._wsave(filename, chain::Chains) = write(filename, chain)

# name2model = Dict(
#     "hospit" => Regional.model_hospit,
#     "hospit2" => Regional.model_hospit_v2,
#     "hospitnp" => Regional.model_hospit_nonparametric,
#     "deaths" => Regional.model_deaths,
#     "cases"  => Regional.model_cases,
#     "cases2"  => Regional.model_cases_v2,
#     "intdeaths" => Regional.model_international_deaths,
# )
#
# model2observable = Dict(
#     "hospit" => ["hospit"],
#     "hospit2" => ["hospit"],
#     "hospitnp" => ["hospit"],
#     "deaths" => ["death"],
#     "cases"  => ["cases"],
#     "cases2"  => ["cases"],
#     "intdeaths" => ["death"],
# )

"Converts a vector of tuples to a tuple of vectors."
function vectup2tupvec(ts::AbstractVector{<:Tuple})
    k = length(first(ts))

    return tuple([[t[i] for t in ts] for i = 1:k]...)
end

"Converts a vector of named tuples to a tuple of vectors."
function vectup2tupvec(ts::AbstractVector{<:NamedTuple})
    ks = keys(first(ts))

    return (; (k => [t[k] for t in ts] for k ∈ ks)...)
end


"""
    rename!(d::Dict, names::Pair...)

Renames the keys given by `names` of `d`.
"""
function rename!(d::Dict, names::Pair...)
    # check that keys are not yet present before updating `d`
    for k_new in values.(names)
        @assert k_new ∉ keys(d) "$(k_new) already in dictionary"
    end

    for (k_old, k_new) in names
        d[k_new] = pop!(d, k_old)
    end
    return d
end

"""
convert a vector (a) of vectors (b) of vectors (c) into a 3 dim Array[a,b,c].
"""
function arrarrarr2arr(a::AbstractVector{<:AbstractVector{<:AbstractVector{T}}}) where {T<:Real}
    n1, n2, n3 = length(a), length(first(a)), length(first(first(a)))

    A = zeros(T, (n1, n2, n3))
    for i = 1:n1
        for j = 1:n2
            for k = 1:n3
                A[i, j, k] = a[i][j][k]
            end
        end
    end

    return A
end

"""
    Chains(d::Dict)

Converts a `Dict` into a `Chains`, assuming the values of `d` is either
- `AbstractVector`: samples for a `Real` variable
- `AbstractMatrix`: samples for a `Vector` variable
- `AbstractArray{<:Any, 3}`: samples for a `Matrix` variable
"""
function MCMCChains.Chains(d::Dict)
    vals = []
    names = []

    for (k, v) in pairs(d)
        if v isa AbstractVector
            push!(vals, v)
            push!(names, k)
        elseif v isa AbstractMatrix
            push!(vals, v)
            append!(names, ["$(k)[$(i)]" for i = 1:size(v, 2)]) # assuming second dimension is dimensionality
        elseif v isa AbstractArray{<:Any, 3}
            indices = CartesianIndices(v[1, :, :])

            # The ordering is such that when you call `reshape` on the vector, the symbols will correspond with
            # the actual indices in the matrix, e.g. `X[i, j]` will be the same as `reshape(X, size)[i, j]`.
            for i = 1:size(indices, 2)
                for j = 1:size(indices, 1)
                    push!(vals, v[:, j, i])
                    push!(names, "$(k)[$j, $i]")
                end
            end
        else
            throw(ArgumentError("I'm so, so sorry but I can't handle $(typeof(v)) :("))
        end
    end

    return Chains(reduce(hcat, vals[2:end]; init = vals[1]), reduce(vcat, names[2:end]; init = names[1]))
end

"""
# hpdi
Compute high density region.
Derived from `hpd` in MCMCChains.jl.
By default alpha=0.11 for a 2-sided tail area of p < 0.055% and p > 0.945%.
"""
function hpdi(x::AbstractVector{T}; alpha=0.11) where {T<:Real}
    n = length(x)
    m = max(1, ceil(Int, alpha * n))

    y = sort(x)
    a = y[1:m]
    b = y[(n - m + 1):n]
    _, i = findmin(b - a)

    return [a[i], b[i]]
end
