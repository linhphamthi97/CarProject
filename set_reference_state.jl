#=
author: Linh Pham Thi
created on: 02.01.2020

Need to use package https://github.com/KristofferC/NearestNeighbors.jl before running.
=#

using Ipopt, JuMP, LinearAlgebra, NearestNeighbors, DelimitedFiles

function set_reference_state(x_current)
    x_ref = Vector{Float64}(undef, 4)

    # Choosing the closest centerline point to current location and use that as reference
        # Load data from csv file, skip first row (header)
        # data_track = readdlm(".//3yp_track2500.csv", ',', Float64, skipstart = 1)

        # Transpose to fit format of NearestNeighbors package and using the first 2 columns only
        # data = transpose(data_track[:,1:2])
        v1 = collect(range(-2, stop=5, length=2500))
        v2= ones(2500)
        data = transpose([v1 v2])

        # Create k-d tree using the data (centerline points)
        tree = KDTree(data, Euclidean(); leafsize=10)

        # Query closest point to current position on centerline
        centerline_point_index, centerline_distance = knn(tree, x_current[1:2], 1, false)
        centerline_point = data[:,centerline_point_index]

        println("Centerline point: ", centerline_point)


    # Set reference state
    x_ref[1] = centerline_point[1]        # x coordinate
    x_ref[2] = centerline_point[2]        # y coordinate
    x_ref[3] = 0                          # angular velocity
    x_ref[4] = 1                          # velocity of centre of mass

    return x_ref

end
