#=
Need to use package https://github.com/KristofferC/NearestNeighbors.jl before running.

This file contains a function to set the reference states given the current state in the (s,d) coordinate system.

=#

using LinearAlgebra, NearestNeighbors, DelimitedFiles
include("settings.jl")
include("track.jl")

function set_reference_state(x_current, x_ref)
    # Choose the centerline point v*dt away from current centerline point
    x_ref[1, 1] = x_current[1] + x_current[4] * settings.delta_t
    x_ref[2, 1] = 0

    for i=2:settings.horizon_length
        x_ref[1, i] = x_ref[1, i-1] + x_current[4] * settings.delta_t
    end

    x_ref[2, :] .= 0                          # centerline deviation
    x_ref[3, :] .= 0                          # angular velocity
    x_ref[4, :] .= settings.target_velocity   # velocity of centre of mass

    return track.data_track
end
