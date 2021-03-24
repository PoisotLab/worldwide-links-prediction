# Prediction of the number of links from the number of species using the flexible links model
# and estimation of the probability of networks being stable according to their species richness

theme(:mute)
default(; frame=:box)

## Fit flexible links model using all food webs archived on Mangal

# Flexible links model and functions
@model FL(S,R) = begin
  N = length(S)
  # Number of trials
  F = S.*S .- (S.-1)
  # Parameters
  ϕ ~ Normal(3.0, 0.5)
  μ ~ Beta(3.0, 7.0)
  for i in 1:N
    R[i] ~ BetaBinomial(F[i], μ*exp(ϕ), (1-μ)*exp(ϕ))
  end
  return μ, ϕ
end

# Get the FL model
ls = DataFrame(CSV.read(joinpath("data", "mangal_foodwebs.csv"), DataFrame))
S = vec(ls[:,:S])
L = vec(ls[:,:L])
R = L .- (S.-1)
chain = sample(FL(S,R), HMC(0.01,10), 3000)

# Beta-binomial for the FL model
μ = mean(get_params(chain).μ)
ϕ = mean(get_params(chain).ϕ)


## Simulate counterfactuals

# Predict the number of links from the number of species
function prediction(chain, S)
  p = get_params(chain[200:end,:,:])
  i = rand(1:length(p.μ))
  μ, ϕ = p.μ[i], p.ϕ[i]
  return rand(BetaBinomial(S^2-(S-1), μ*exp(ϕ), (1-μ)*exp(ϕ))) + (S-1)
end

# Simulate the number of links for a range of species richness
sp = collect(5:1000) # numbers of species
nsp = length(sp)
nsim = 1000 # number of simulations

predicted_links = DataFrame()

for i in 1:nsim
  colname = "Lhat$i"
  Lhat = [prediction(chain, S) for S in sp]
  predicted_links[!, colname] = Lhat
end

## Summary statistics for the predicted numbers of links

measures = DataFrame()

# Species richnes
measures.S = sp

# Median predicted numbers of links
measures.L_median = median.(eachrow(predicted_links))

# 68% percentile intervals of predicted numbers of links
lower_links = quantile.(eachrow(predicted_links), 0.16)
upper_links = quantile.(eachrow(predicted_links), 0.84)
measures.L_PI = upper_links .- lower_links


## Estimation of the probability of the networks being stable (stability score)

# Integration over σ and get the area under the curve
function ∫(x,y)
  @assert length(x) == length(y)
  S = zero(Float64)
  for i in 2:length(x)
    S += (x[i]-x[i-1])*(y[i]+y[i-1])*0.5
  end
  return S
end

# This function is basically a wrapper around the code before
function stability_score(S, μ, ϕ)
    # This is the flexible link distribution
    D = BetaBinomial(S*S-(S-1), μ*exp(ϕ), (1-μ)*exp(ϕ))
    # We know that 0 < σ ≤ √(S/(S-1))
    Σ = LinRange(1e-5, sqrt(S/(S-1)), 100)
    # This is to store the values of maximal links before the network gets unstable
    V = zeros(Float64, length(Σ))
    for (i,σ) in enumerate(Σ)
        # This is the number of links for which the network is unstable at a given S and σ
        threshold = S/(σ*σ)-(S-1)
        # We can get the probability of the network being unstable with the CDF
        V[i] = cdf(D, threshold)
    end
    return ∫(Σ, V)
end

# Stability of networks of 5 to 1000 species
stab = zeros(Float64, length(sp))
p = Progress(length(sp))
Threads.@threads for i in 1:length(sp)
    stab[i] = stability_score(sp[i], μ, ϕ)
    next!(p)
end

measures.stab = stab

# Stability score as a function of species richness
plot(sp, stab, lab="", fill=(0, 0.1))
xaxis!(:log, extrema(sp), "Species richness")
yaxis!((0, 1), "Stability score")
savefig(joinpath("figures", "stability_score.png"))

# Write file
CSV.write(joinpath("data", "measures.csv"), measures)
