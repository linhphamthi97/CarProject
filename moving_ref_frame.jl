#=
author: Linh Pham Thi
created on: 30.03.2020
=#
include("settings.jl")
include("track_s.jl")
using Plots, DelimitedFiles, NearestNeighbors

# =============================================================================
# Track module
# =============================================================================
module track_s
export sections, datapoints, section_number
include("settings.jl")
include("track_s.jl")
using NearestNeighbors
#=
This module serves as the base data for the track information. It specifies the
datapoints and degrees based on which the B-spline curves will be evaluted.
Two dictionaries are part of this module:
datapoints:
    - number of entries is either 1 (substitute track) or 17 (original)
    - key is the section number
    - value: datapoints array, degree for B-spline

sections:
    - number of entries based on datapoints
    - key is the section number
    - values: datapoints on the centerline in (x, y) coordinates) (B-spline curve evaluated)
              a tree containing all of the datapoints on the centerline in the section
              L - length of the section
              P - control points of the section
              U - knot vector of the section
              [s,d] - corresponding (s,d) coordinates
              a tree containing all of the centerline points (s,d) coordinates
=#
# ..............................................................................
# Datapoints and degrees (degree of B-spline curve)
    datapoints = Dict{Int64, Any}()
    if settings.original_track == true
        # Original track
        datapoints[1] = ([0 1; 2.391 1], 1)
        datapoints[2] = ([2.391 1; 2.805 0.914; 3 0.815; 3.403 0.732; 3.879 0.854; 4.033 0.93; 4.405 1], 3)
        datapoints[3] = ([4.405 1; 6.789 1], 1)
        datapoints[4] = ([6.789 1; 7.225 0.823; 7.398 0.444; 7.222 -0.026; 7 -0.166], 2)
        datapoints[5] = ([7 -0.166; 4.823 -0.616], 1)
        datapoints[6] = ([4.823 -0.616; 4.2 -1; 4 -1.543], 2)
        datapoints[7] = ([4 -1.543; 3.794 -3.467], 1)
        datapoints[8] = ([3.794 -3.467; 3.678 -3.687; 3.385 -3.8; 3.117 -3.682; 3.009 -3.483], 2)
        datapoints[9] = ([3.009 -3.483; 2.6 -0.575], 1)
        datapoints[10] = ([2.6 -0.575; 2.446 -0.405; 2.226 -0.461], 2)
        datapoints[11] = ([2.226 -0.461; 1.507 -1.258], 1)
        datapoints[12] = ([1.507 -1.258; 1.427 -1.423; 1.4 -1.605], 2)
        datapoints[13] = ([1.4 -1.605; 1.4 -3.033], 1)
        datapoints[14] = ([1.4 -3.033; 1.179 -3.551; 0.649 -3.8; 0.0246 -3.555; -0.172 -3.209], 2)
        datapoints[15] = ([-0.172 -3.209; -0.2 -3.015; -0.339 -2.549; -0.619 -2.296; -0.839 -2.096; -1 -1.589], 3)
        datapoints[16] = ([-1 -1.589; -1 0], 1)
        datapoints[17] = ([-1 0 ; -0.64 0.768; 0.0 1.0], 2)

    else # Substitute track
        datapoints[1] = ([0 0; 1 1.5 ; 2.5 2.0 ; 4 1.5; 5.5 1; 7 2 ; 8 2.5; 9.5 0;
              8 -2; 7 -1.5; 5.5 0 ; 4 -1.5 ; 2.5 -2; 1 -1.5; 0 0], 2)
    end

# .............................................................................
# Create sections, section lengths and binary search tree for each section using data
    sections = Dict{Int64, Any}()
    for i in 1:size(length.(collect(keys(datapoints))),1)
        # Section
        section_data, P, U = B_spline(datapoints[i][1], datapoints[i][2], settings.sample_number)

        # Total section length = integral of derivative from t=0 to 1
        L = length_calculator(P, size(datapoints[i][1],1)-1, datapoints[i][2], U, 0, 1)

        # s-value calculator
        s = zeros(settings.sample_number, 1)
            # Previous sections lengths added
            for k in 1:size(length.(collect(keys(datapoints))),1)-1
                s .+= sections[k][3]
            end

            # Current section up to t
            n = 1
            for t in LinRange(0,1,settings.sample_number)
                s[n] +=  length_calculator(P, size(datapoints[i][1], 1) - 1, datapoints[i][2], U, 0, t)
                n += 1
            end

        sd_coordinates = [s zeros(settings.sample_number, 1)]

        # Add components to dictionary entry
        sections[i] = (section_data, KDTree(transpose(section_data), Euclidean(); leafsize=10),
                    L, P, U, sd_coordinates, KDTree(transpose(sd_coordinates), Euclidean(); leafsize=10))
    end

# .............................................................................
# Global parameters
    section_number = size(length.(collect(keys(track_s.datapoints))),1)

    global total_lap_length = 0
    for i in 1:section_number
        global total_lap_length += sections[i][3]
    end
end

# =============================================================================
# Functions
# =============================================================================
function find_nearest_centerline_point(P_xy)
    distances = Array{Float64}(undef, track_s.section_number, 1)
    t_index = Array{Int64}(undef, track_s.section_number, 1)

    for i in 1:track_s.section_number
        t_value, distance = knn(track_s.sections[i][2], P_xy, 1, false)
        t_index[i] = t_value[1]
        distances[i] = distance[1]
    end

    d, min_index = findmin(distances)
    t_vector = LinRange(0,1, settings.sample_number)
    t = t_vector[t_index[min_index]]

    return t, d, min_index[1]      # curve parameter t, distance d, section number (min_index)
end

function transform_coordinate(P_xy)
    # This function transforms a point P(x,y) into the moving coordinate system P(s,d)
    # Input: P - point with x and y coordinate
    # Output: s, d coordinates of point

    # Find nearest centerline point
    t, d, section_index = find_nearest_centerline_point(P_xy)

    # Compute s by adding together section lengths and partial section length
        s = 0

        # Previous section lengths
        for i in 1:(section_index-1)
            s += track_s.sections[i][3]
        end

        # Partial section length
        s += length_calculator(track_s.sections[section_index][4], size(track_s.datapoints[section_index][1],1)-1, track_s.datapoints[section_index][2], track_s.sections[section_index][5], 0, t)

    # determine sign of d (left - positive or right - negative of centerline)
    derivatives = derivatives_calculator(track_s.sections[section_index][4], size(track_s.datapoints[section_index][1],1)-1,
                                         track_s.datapoints[section_index][2], track_s.sections[section_index][5], t)

        # d positive
        xy_centerline = transform_coordinate_back([s;0])
        x1 = xy_centerline[1] - d*derivatives[2]
        y1 = xy_centerline[2] + d*derivatives[1]
        dist1 = (x1 - P_xy[1])^2 + (y1 + P_xy[2])^2

        # d negative
        x2 = xy_centerline[1] + d*derivatives[2]
        y2 = xy_centerline[2] - d*derivatives[1]
        dist2 = (x2 - P_xy[1])^2 + (y2 + P_xy[2])^2

        if dist2<dist1
            d = -d
        end

    return s, d
end

function transform_coordinate_back(P_sd)
    # This function transforms a point P(s,d) into the fixed coordinate system P(x,y)
    # Input: s, d coordinates of point
    # Output: P - point with x and y coordinate

    # Multiple laps case
    while P_sd[2] >= track_s.total_lap_length
        P_sd[2] -= track_s.total_lap_length
    end

    # Find nearest centerline point and it's x-y coordinate
        tree_distances = Array{Float64}(undef, track_s.section_number, 1)
        sd_index = Array{Int64}(undef, track_s.section_number, 1)

        for i in 1:track_s.section_number
            sd_value, distance = knn(track_s.sections[i][7], P_sd, 1, false)
            sd_index[i] = sd_value[1]
            tree_distances[i] = distance[1]
        end

        d_min, min_idx = findmin(tree_distances)

        min_index = min_idx[1]
                    # Selecting relevant section, 1 is the (x,y) coordinate bit, [sd_index[min_index]] is the index to the relevant point
        xy_centerline = track_s.sections[min_index][1][sd_index[min_index], :]

    t = LinRange(0,1, settings.sample_number)[sd_index[min_index]]

    # Use derivative and d to add deviation bit
    derivatives = derivatives_calculator(track_s.sections[min_index][4], size(track_s.datapoints[min_index][1],1)-1,
                                         track_s.datapoints[min_index][2], track_s.sections[min_index][5], t)
    norm = sqrt(derivatives[1]^2 + derivatives[2]^2)
    x = xy_centerline[1] - P_sd[2]*derivatives[2]/norm
    y = xy_centerline[2] + P_sd[2]*derivatives[1]/norm

    return [x; y]
end

function curvature(s,d)     # Input could be done as a vector, but needs to be separate for the JuMP user defined function
    # Find nearest centerline point and it's section data
        tree_distances = Array{Float64}(undef, track_s.section_number, 1)
        sd_index = Array{Int64}(undef, track_s.section_number, 1)

        for i in 1:track_s.section_number
            sd_value, distance = knn(track_s.sections[i][7], [s;d], 1, false)
            sd_index[i] = sd_value[1]
            tree_distances[i] = distance[1]
        end

        d_min, min_idx = findmin(tree_distances)

        section_index = min_idx[1]

        t = LinRange(0,1, settings.sample_number)[sd_index[section_index]]

    # Find curvature
        derivatives = derivatives_calculator(track_s.sections[section_index][4],
                                       size(track_s.datapoints[section_index][1],1) - 1,
                                       track_s.datapoints[section_index][2],
                                       track_s.sections[section_index][5],
                                       t)
    return derivatives[5]
end
# =============================================================================
# Plot generation
# =============================================================================
#=
if settings.original_track == true
    data_track = Array{Float64}(undef, 2500,3)
    data_track[:,1:3] = readdlm(".//3yp_track2500.csv", ',', Float64, skipstart = 1)

    p = plot(data_track[:,1], data_track[:,2], color=:black, xlim=(-2, 8), ylim=(-4, 2), label="Data", legend=:bottomright)
    # p = plot(data_track[:,1], data_track[:,2], color=:black, xlim=(-1.2, 0.5), ylim=(-4, -1), label="", aspect_ratio=:equal, legend=:bottomright)
else
    p = plot(xlim=(-0.5, 10), ylim=(-2.5, 3), label="", aspect_ratio=:equal)
end

for i in 1:size(length.(collect(keys(track_s.datapoints))),1)
    # Curve
    plot!(p, track_s.sections[i][1][:,1], track_s.sections[i][1][:,2], label="")

    # Datapoints
    plot!(p, track_s.datapoints[i][1][:,1], track_s.datapoints[i][1][:,2], color=:purple, seriestype = :scatter, label="")
end
display(p)
savefig("e:\\4yp_images\\figure1")

A = transform_coordinate_back([0; 0.03])
println("A= ", A)

B = transform_coordinate(A)

println("B= ", B)

=#
