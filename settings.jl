#=
author: Linh Pham Thi
created on: 02.02.2020

This file contains all of the settings to change about the simulation or the model.
=#

module settings
export u_initial, x_initial, horizon_length, delta_t, sample_number,
        model_type, original_track, track_width, N_sim, target_velocity,
        Q_v, R_laptime
using LinearAlgebra
# Track selector & parameters
    original_track = false      # true for original track, false for substitute track
    track_width = 0.18          # Track width (total) for alternative track

# Simulation parameters
    delta_t = 0.01
    model_type = 1      # Bicycle model type to use in the simulation
                        # 1 for kinematic, 2 for dynamic (doesn't work)
    N_sim = 30000         # number of iterations in the simulation
    target_velocity = 3

# Initial conditions
    # (x,y) coordinate system
    # u_initial = [0., 0.]
    # x_initial = [0, 1., 0.1, 1]

    # (s,d) coordinate system
    u_initial = [0., 0.]
    x_initial = [0, 0.0, 0, 0.001]

# Control model (bicycle model) parameters
    horizon_length = 20

    # Reference state tracker


    # Laptime minimiser
    Q_v = 1000
    R_laptime = Diagonal(vec(fill(4,2)))

# (s,d) coordinate system
    sample_number = 1000        # Number of samples (different t) on a section
    integral_steps = 1000       # 1/dt = integral_steps for length calculation

end # module
