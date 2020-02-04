#=
author: Linh Pham Thi
created on: 02.01.2020

This is the model for the car's behaviour for the simulation. This model takes
a current state and the inputs applied to output the next state. The model might
be the same as the one used in the MPC but ideally it should be a more elaborate
model to truly simulate the discrepancies between the simplified model used in
the controller and the real world.
=#

function simulation(x, u)
    throttle = u[1,1]
    delta = u[2,1]
    phi = x[3,1]
    v = 1 #x[4,1]

    delta_t = 0.1

    x_next = Vector{Float64}(undef, 4)

    # Bicycle model and track boundary constraints
    # Steering slip - negative for oversteer, positive for understeer
    C1 = 0.5
    # Steering angle coupling - radians turned / meter travelled
    C2 = 15
    # Duty cycle to acceleration a = Cm1 * throttle
    Cm1 = 19.6
    # Cm2 = Cm1 / v_motor_max (i.e. max speed with no air resistance)
    Cm2 = 14.2
    # Rolling resistance
    Cr1 = 1.8
    # Reduced air resistance coefficient (0.5 * rho * A * C_d)
    Cr2 = 0.0045

# term for v_dot: (Cm1 * throttle - Cm2 * throttle * v - Cr2 * v * v - Cr1 * sign(v)))

    x_next = x + [v * cos(phi + C1 * delta) * delta_t;
                 v * sin(phi + C1 * delta)* delta_t;
                 C2 * delta * v * delta_t;
                 0]

    return x_next
end
