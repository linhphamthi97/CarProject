#=
author: Linh Pham Thi
created on: 02.01.2020

This is the main file of the simulation.
The structure of the program is as follows:
1. Calculate a target trajectory / target states 'offline'
2. while simulation is running:
    2.1. Check current state and generate state reference
    2.2. Apply MPC, subject to model contraints (bicycle model or combined slip model)
    2.3. Apply the input suggested by the MPC and simulate the next state
3. Plot results
=#

using Plots, Ipopt, JuMP, DelimitedFiles, CSV
include("run_mpc.jl")
include("build_mpc.jl")
include("set_reference_state.jl")
include("simulation.jl")
include("settings.jl")

function main()
    # Use the GR backend for Plots.jl, because it's fast
    gr()

    # Intialise parameters

    N_sim = 50
    u = settings.u_initial
    x = settings.x_initial
    x_ref = Array{Float64}(undef, 4, settings.horizon_length)
    x_plan = repeat(settings.x_initial, 1, settings.horizon_length)
    u_plan = repeat(settings.u_initial, 1, settings.horizon_length)


    # just keep the history for debugging
    x_hist = zeros(Float64, N_sim + 1, 4)
    x_hist[1, :] = x
    u_hist = zeros(Float64, N_sim, 2)

    # initialise the MPC solver
    mpcsolver = build_mpc(settings.horizon_length)


    # Create an animation
    anim = @animate for i in 1:N_sim

        # Plot the current position and the track
        plot([x[1]], [x[2]], marker=(:hex, 10), xlim=(-2, 8), ylim=(-4, 2))

        # Set reference state
        track_data = set_reference_state(x, x_ref)

        # Run the MPC control optimization
        x_plan, u_plan = run_mpc(mpcsolver, x, x_ref, x_plan, u_plan, settings.horizon_length)

        # extract the first input
        u .= u_plan[:,1]

        # Draw the planned future states from the MPC optimization
        plot!(x_plan[1, :], x_plan[2, :], linewidth=5, label="Predicted state")
        plot!(track_data[1,:], track_data[2,:], label="Track centerline")

        # Apply the planned inputs and simulate one step in time
        x .= simulation(x, u)

        # save the history of states and inputs
        x_hist[i + 1, :] = x
        u_hist[i, :] = u

        # println("x after simulation= ", x)
        # println(" ")
        # println(" ")
    end

    # Save animation as a gif file
     #gif(anim, "/mpc.gif", fps = 20)

    return x_plan, x_ref, u_plan, x_hist, u_hist
end

x_plan, x_ref, u_plan, x_hist, u_hist = main()

# plot the race track and states
track2500 = CSV.read("3yp_track2500.csv");
plot(track2500.x, track2500.y, label = "racetrack-centerline")
plot!(x_hist[:, 1], x_hist[:, 2], label = "car trajectory")
