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

using Plots, Ipopt, JuMP
include("run_mpc.jl")
include("set_reference_state.jl")
include("simulation.jl")
include("settings.jl")

function main()
    # Use the GR backend for Plots.jl, because it's fast
    gr()

    # Intialise parameters
    u = settings.u_initial
    x = settings.x_initial
    x_ref = Vector{Float64}(undef, 4)
    x_plan = repeat(settings.x_initial, 1, 2)
    u_plan = Vector{Float64}(undef, 2)

    # Create an animation
    anim = @animate for i in 1:1
        println("Timestep ", i)

        # Plot the current position and the track
        plot([x[1]], [x[2]], marker=(:hex, 10), xlim=(-2, 8), ylim=(-4, 2))

        # Set reference state
        x_ref = set_reference_state(x)
        # println("Current state: ", x)
        # println(" Reference state: ", x_ref)

        # Run the MPC control optimization
        x_plan, u = run_mpc(x_ref, x, u, x_plan, settings.horizon_length)

        # Draw the planned future states from the MPC optimization
        plot!(x_plan[1, :], x_plan[2, :], linewidth=5)

        # Apply the planned inputs and simulate one step in time
        x = simulation(x, u)
        # println("x after simulation= ", x)
        # println(" ")
        # println(" ")

    end

    # Save animation as a gif file
    gif(anim, "/mpc.gif", fps = 10)

end

main()
