using JLD2, MAT, MCMCChains, RemoteFiles

# Download results from https://github.com/ap-browning/internalisation
jld2_file = @RemoteFile(Data,
    "https://github.com/ap-browning/internalisation/raw/main/Results/MainResult.jld2",
    dir = @__DIR__, file = "Browning2021.jld2")
download(jld2_file)

# Load results
@load "Browning2021.jld2"

# Convert chains to a matrix
Base.Matrix(C::Chains) = permutedims(hcat([Matrix(C[:,i,:])[:] for i = 1:size(C,2)-1]...))
C_Matrix = collect(Matrix(C)')

# Save to a .mat file
mat_file = matopen("Browning2021.mat", "w")
write(mat_file, "theta", C_Matrix)
write(mat_file, "theta_best", θ)
write(mat_file, "param_names", String.(C_pilot.name_map.parameters))
close(mat_file)

