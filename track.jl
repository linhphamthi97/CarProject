#=
author: Linh Pham Thi
created on: 24.03.2020

This file contains srcipt to generate the track data needed (centerline, inner
and outer track bounds, converting to (s,a) based coordinate system).

data_track array columns:
1, 2: (x, y) coordinate of centerline
3: width of track
4, 5: (x, y) coordinate of inner track bound
6, 7: (x, y) coordinate of outer track bound

=#
include("settings.jl")

module track
export data_track, data_centerline, tree
using LinearAlgebra, NearestNeighbors, DelimitedFiles

    # Load data from csv file, skip first row (header)
    data_track = Array{Float64}(undef, 2500,7)
    data_track[:,1:3] = readdlm(".//3yp_track2500.csv", ',', Float64, skipstart = 1)

    # Creating the track bound coordinates
    for i in 1:size(data_track, 1)
        if i == 1
            vector1 = data_track[size(data_track, 1)-1, 1:2]
            vector2 = data_track[i+1, 1:2]
        elseif i == size(data_track, 1)
            vector1 = data_track[i-1, 1:2]
            vector2 = data_track[1, 1:2]
        else
            vector1 = data_track[i-1, 1:2]
            vector2 = data_track[i+1, 1:2]
        end

        tangent = normalize(vector2 - vector1)

        orthogonal1 = [- tangent[2]; tangent[1]]
        orthogonal2 = - orthogonal1

        # inner bound
        data_track[i,4:5] = -orthogonal1 * data_track[i,3]/2 + data_track[i,1:2]

        # outer bound
        data_track[i,6:7] = orthogonal1 * data_track[i,3]/2 + data_track[i,1:2]
    end

    # TODO: transformation to (s,a) coordinate system

    # Transpose to fit format of NearestNeighbors package and using the first 2 columns only
    data_centerline = transpose(data_track[:,1:2])

    # Create k-d tree using the data (centerline points)
    tree = KDTree(data_centerline, Euclidean(); leafsize=10)

end
