#=
author: Linh Pham Thi
created on: 02.02.2020

This file contains all of the settings to change about the simulation or the model.
=#

module settings
export u_initial, x_initial, horizon_length, delta_t
# Simulation parameters
    delta_t = 0.1

# Initial conditions
    u_initial = [0, 0]
    x_initial = [0, 1, 0.1, 1]

# Control model (bicycle model) parameters
    horizon_length = 10



end # module
