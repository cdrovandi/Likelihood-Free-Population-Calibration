## INTERNALISATION MODEL POSTERIOR PREDICTIVE
#
#  Calls code available at https://github.com/ap-browning/internalisation (run from the root directory)

using Inference
using Model

using JLD2
#using MAT
using KernelDensity
using Plots
using StatsFuns
using StatsPlots
using Random

# Load model and results
include("Results/MainResult_Setup.jl")
@load "Results/MainResult.jld2"

# 1000 posterior samples
Θ = sample(C[1:100:end],1000)
x = range(0.0,1.0,length=200)

# Some observation times...
t = [1.0,5.0,30.0,180.0]

# Store h ← pdf of fraction internalised as a function of (x,t)
h1 = zeros(length(x),length(Θ))
h2 = zeros(length(x),length(Θ))
h3 = zeros(length(x),length(Θ))
h4 = zeros(length(x),length(Θ))
H = [h1,h2,h3,h4]
B = zeros(length(x),4)
for (i,tᵢ) in enumerate(t)
    # for (j,θⱼ) in enumerate(Θ)
    #     d = param_dist(θⱼ[5:end-1])
    #     ξ = sample(d,100_000)
    #     Z = parameters_to_antibody(tᵢ,ξ[2,:],ξ[3,:],θ[end])
    #     H[i][:,j] = pdf(kde(Z[2] ./ Z[1]),x)
    # end
    # At best-fit
    ξ = sample(param_dist(θ[5:end-1]),100_000)
    Z = parameters_to_antibody(tᵢ,ξ[2,:],ξ[3,:],θ[end])
    B[:,i] = pdf(kde(Z[2] ./ Z[1]),x)
end


# mat_file = matopen("Browning2021_Predictive.mat", "w")
# write(mat_file, "x", x)
# write(mat_file, "t", t)
# write(mat_file, "h1", h1)
# write(mat_file, "h2", h2)
# write(mat_file, "h3", h3)
# write(mat_file, "h4", h4)
# write(mat_file, "B", B)
# close(mat_file)