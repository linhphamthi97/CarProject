#=
author: Linh Pham Thi
created on: 02.01.2020

Need to use package https://github.com/KristofferC/NearestNeighbors.jl before running.
=#

using Ipopt, JuMP, LinearAlgebra, NearestNeighbors, DelimitedFiles
include("settings.jl")

function set_reference_state(x_current)
    x_ref = Array{Float64}(undef, 4, settings.horizon_length)

    # Choosing the closest centerline point to current location and use that as reference
        # Load data from csv file, skip first row (header)
        data_track = readdlm(".//3yp_track2500.csv", ',', Float64, skipstart = 1)

        # Transpose to fit format of NearestNeighbors package and using the first 2 columns only
        data = transpose(data_track[:,1:2])

        # Create k-d tree using the data (centerline points)
        tree = KDTree(data, Euclidean(); leafsize=10)

        # Query closest point to current position on centerline
        centerline_point_index, centerline_distance = knn(tree, x_current[1:2], 1, false)
        centerline_point = data[:,centerline_point_index]

        # Set reference state
        x_ref[1,1] = centerline_point[1]        # x coordinate
        x_ref[2,1] = centerline_point[2]        # y coordinate

        # Query next points, within distance v*delta_t away
        for i=2:settings.horizon_length
            idx = inrange(tree, x_ref[1:2,i-1], x_current[4] * settings.delta_t, false)
            x_ref[1:2,i] = data[:, maximum(idx)]
        end

    # Set reference state
    x_ref[3, :] .= 0                          # angular velocity
    x_ref[4, :] .= 1                          # velocity of centre of mass

    return x_ref

end
