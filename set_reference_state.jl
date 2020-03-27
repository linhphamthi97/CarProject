#=
author: Linh Pham Thi
created on: 02.01.2020

Need to use package https://github.com/KristofferC/NearestNeighbors.jl before running.

This file contains a function to set the reference states given the current state.

=#

using LinearAlgebra, NearestNeighbors, DelimitedFiles
include("settings.jl")
include("track.jl")

function set_reference_state(x_current, x_ref)
    # Choosing the closest centerline point to current location and use that as reference

        # Query closest point to current position on centerline
        close_point_index, distances = knn(track.tree, x_current[1:2], 2, false)
        close_points = track.data_centerline[:,close_point_index]

        # Linear interpolation between 2 closest points to get centerline point
        r = sum((close_points[1,i] - close_points[2,i])^2 for i in 1:2)
        r1 = (distances[2]^2 - distances[1]^2 + r)/(2*r)
        r2 = 1- r1

        # Set reference state
        x_ref[1,1] = close_points[1,1]*r1 + close_points[1,2]*r2        # x coordinate
        x_ref[2,1] = close_points[2,1]*r1 + close_points[2,2]*r2        # y coordinate

        # Query next points, within distance v*delta_t away
        for i=2:settings.horizon_length
            idx = inrange(track.tree, x_ref[1:2,i-1], x_current[4] * settings.delta_t, true)

            # If points within that range are equi-spaced, then choose the one
            # with the highest index, i.e. the point of distance v*delta_t away
            # in the direction of the car moving (clockwise)
            # If there is a discontinuity, then
            #       if the discontinuity is big (i.e. point where the track data
            #           begins) or if discontinuity is after the index of
            #           x_ref[1:2,i-1] point: use the point before discontinuity
            #       else use the point with highest index

            if maximum(idx[2:end]-idx[1:end-1]) - minimum(idx[2:end]-idx[1:end-1]) == 0
                data_idx = maximum(idx)
            else
                discont_value, discont_idx = findmax(idx[2:end]-idx[1:end-1])

                if discont_value>=0.9*size(track.data_centerline, 2) || idx[discont_idx] > data_idx
                    data_idx = idx[discont_idx]
                else
                    data_idx= maximum(idx)
                end

            end

            x_ref[1:2,i] = track.data_centerline[:, data_idx]

        end

    # Set reference state
    x_ref[3, :] .= 0                          # angular velocity
    x_ref[4, :] .= 1                          # velocity of centre of mass

    return track.data_track
end
