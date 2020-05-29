#=
author: Linh Pham Thi

This is the main file of the simulation for the models in the (s,d) coordinate
system.
=#

using Plots, Ipopt, JuMP, DelimitedFiles, CSV
include("simulation_sd.jl")
include("run_mpc.jl")
include("build_mpc_sd_reference_tracker.jl")
include("set_reference_state_sd.jl")
include("settings.jl")
include("track_s.jl")
include("moving_ref_frame.jl")

function main()

    # Use the GR backend for Plots.jl, because it's fast
    gr()

    # Intialise parameters
    u = settings.u_initial
    x = settings.x_initial
    x_ref = Array{Float64}(undef, 4, settings.horizon_length)
    x_plan = repeat(settings.x_initial, 1, settings.horizon_length)
    x_plan_xy = Array{Float64}(undef, 4, settings.horizon_length)
    u_plan = repeat(settings.u_initial, 1, settings.horizon_length)
    track_data = Array{Float64}(undef, 2500, 7)
    P_xy = zeros(2,1)


    # just keep the history for debugging
    x_hist = zeros(Float64, settings.N_sim + 1, 4)
    x_hist[1, :] = x
    u_hist = zeros(Float64, settings.N_sim, 2)

    # initialise the MPC solver
    mpcsolver = build_mpc(settings.horizon_length)

    # Create an animation
    # anim = @animate
    for i in 1:settings.N_sim
        println("")
        println("Iteration ", i)
        println("Current state: ", x[:,1])

        # Plot the current position and the track
        # P_xy = transform_coordinate_back(x[1:2])
        # plot(P_xy[1], P_xy[2], marker=(:hex, 10), xlim=(-0.5, 10), ylim=(-2.5, 3), label="", aspect_ratio=:equal)

        # Set reference state
        track_data = set_reference_state(x, x_ref)
        # println("Reference state: ", x_ref)

        # Run the MPC control optimization
        x_plan, u_plan = run_mpc(mpcsolver, x, x_ref, x_plan, u_plan, settings.horizon_length)
        # println("planned x= " , x_plan)
        # println("plan - reference= ", x_plan-x_ref)
        println("first planned u= ", u_plan[:,1])
        # println("First obj sum: ", sum(Q[j,j]*(x_plan[j,i]-x_ref[j,i])^2 for i in 1:settings.horizon_length for j in 1:4))

        # extract the first input
        u .= u_plan[:,1]

        # Draw the planned future states from the MPC optimization
        # for i in 1:settings.horizon_length
        #     x_plan_xy[i,:] = transform_coordinate_back(x_plan[i,1:2])
        # end
        # plot!(x_plan_xy[1, :], x_plan_xy[2, :], linewidth=5, label="Predicted state")
        # plot!(track_data[:,1], track_data[:,2], label="Track centerline", color=:blue)
        # plot!(track_data[:,4], track_data[:,5], color=:black, label="")
        # plot!(track_data[:,6], track_data[:,7], color=:black, label="")


        # Apply the planned inputs and simulate one step in time
        x .= simulation(x, u)

        # save the history of states and inputs
        x_hist[i + 1, :] = x
        u_hist[i, :] = u

        # Stop simulation just before it reaches the end point
        if x[1] >= 0.99*track_s.total_lap_length
            return x_plan, x_ref, u_plan, x_hist, u_hist, track_data, i*settings.delta_t
        end

    end

    # Save animation as a gif file
    # gif(anim, "/mpc_sd.gif", fps = 20)

    return x_plan, x_ref, u_plan, x_hist, u_hist, track_data, settings.N_sim*settings.delta_t
end

x_plan, x_ref, u_plan, x_hist, u_hist, track_data, laptime = main()
println("Laptime was ", laptime, "s")

# =============================================================================
# Plots
# =============================================================================
# Centerline
    # track2500 = CSV.read("3yp_track2500.csv");            # Original track
    p = plot(track_s.sections[1][1][:,1], track_s.sections[1][1][:,2],
             label = "racetrack-centerline", color=:blue, xlim=(-0.5, 10),
             ylim=(-2.5, 3.5), aspect_ratio=:equal, legend=:topleft)

# Car trajectory
    x_hist_xy=Array{Float64}(undef, settings.N_sim, 2)
    for i in 1:settings.N_sim
        x_hist_xy[i,:] = transform_coordinate_back(x_hist[i,1:2])
    end
    plot!(p, x_hist_xy[:, 1], x_hist_xy[:, 2], label = "car trajectory", color=:red)

# Track boundaries
    boundaries = Array{Float64}(undef, settings.sample_number, 4)
    for i in 1:settings.sample_number
        if i == 1
            vector1 = track_s.sections[1][1][end, :]
            vector2 = track_s.sections[1][1][i+1, :]
        elseif i == settings.sample_number
            vector1 = track_s.sections[1][1][i-1, :]
            vector2 = track_s.sections[1][1][1, :]
        else
            vector1 = track_s.sections[1][1][i-1, :]
            vector2 = track_s.sections[1][1][i+1, :]
        end

        tangent = normalize(vector2 - vector1)

        orthogonal1 = [- tangent[2]; tangent[1]]
        orthogonal2 = - orthogonal1

        # inner bound
        boundaries[i ,1:2] = -orthogonal1 * settings.track_width/2 + track_s.sections[1][1][i, :]

        # outer bound
        boundaries[i,3:4] = orthogonal1 * settings.track_width/2 + track_s.sections[1][1][i, :]
    end

    plot!(p, boundaries[:,1], boundaries[:,2], color=:black, label="")
    plot!(p, boundaries[:,3], boundaries[:,4], color=:black, label="")

savefig("e:\\4yp_images\\figure")
display(p)

# Curvatures plot
curvatures = Array{Float64}(undef, settings.sample_number,1)
for i in 1:settings.sample_number
    curvatures[i] = curvature(track_s.sections[1][6][i,1], 0)
end

p1 = plot(track_s.sections[1][6][1:999,1], curvatures[1:999], label="", xlabel="Centerline distance s", ylabel="Curvature")

savefig("e:\\4yp_images\\curvature")
