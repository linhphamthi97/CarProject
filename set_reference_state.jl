#=
author: Linh Pham Thi
created on: 02.01.2020

Need to use package https://github.com/KristofferC/NearestNeighbors.jl before running.
=#

using Ipopt, JuMP, LinearAlgebra, NearestNeighbors, DelimitedFiles
include("settings.jl")

# Choosing the closest centerline point to current location and use that as reference
# Load data from csv file, skip first row (header)
data_track = readdlm(".//3yp_track2500.csv", ',', Float64, skipstart = 1)
# Transpose to fit format of NearestNeighbors package and using the first 2 columns only
data = transpose(data_track[:,1:2])

# Create k-d tree using the data (centerline points)
tree = KDTree(data, Euclidean(); leafsize=10)

function set_reference_state(x_current, x_ref)

        # Query closest point to current position on centerline
        centerline_point_index, centerline_distance = knn(tree, x_current[1:2], 1, false)
        centerline_point = data[:,centerline_point_index]

        # Set reference state
        x_ref[1,1] = centerline_point[1]        # x coordinate
        x_ref[2,1] = centerline_point[2]        # y coordinate

        data_idx=centerline_point_index[1]

        # Query next points, within distance v*delta_t away
        for i=2:settings.horizon_length
            idx = inrange(tree, x_ref[1:2,i-1], x_current[4] * settings.delta_t, true)

            # If points within that range are equi-spaced, then choose the one
            # with the highest index, i.e. the point of distance v*delta_t away
            # in the direction of the car moving (clockwise)
            # If there is a discontinuity, then
            #       if the discontinuity is big (i.e. point where the track data
            #           begins) or if discontinuity is after the index of
            #           x_ref[1:2,i-1] point: use the point before discontinuity
            #       else use the point with highest index

            if maximum(idx[2:end]-idx[1:end-1]) - minimum(idx[2:end]-idx[1:end-1]) == 0
                data_idx= maximum(idx)
            else
                discont_value, discont_idx = findmax(idx[2:end]-idx[1:end-1])

                if discont_value>=0.9*size(data,2) || idx[discont_idx] > data_idx
                    data_idx = idx[discont_idx]
                else
                    data_idx= maximum(idx)
                end

            end

            x_ref[1:2,i] = data[:, data_idx]

        end

    # Set reference state
    x_ref[3, :] .= 0                          # angular velocity
    x_ref[4, :] .= 1                          # velocity of centre of mass

    return data
end
