#=
author: Linh Pham Thi
created on: 02.01.2020

This is the model for the car's behaviour for the simulation. This model takes
a current state and the inputs applied to output the next state. The model might
be the same as the one used in the MPC but ideally it should be a more elaborate
model to truly simulate the discrepancies between the simplified model used in
the controller and the real world.

This is the version for the (s,d) coordinate frame.
=#

function simulation(x, u)
    throttle = u[1,1]
    omega = u[2,1]
    e_phi = x[3,1]
    v = x[4,1]
    kappa = curvature(x[1,1], x[2,1])

    x_next = Vector{Float64}(undef, 4)

    # Parameters
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

    x_next = x + [v * cos(e_phi + C1 * omega)/(1-kappa*x[2, 1]);
        v * sin(e_phi + C1 * omega);
        C2 * omega - (kappa * v * cos(e_phi + C1 * omega)/(1-kappa * x[2, 1]));
        (Cm1 - v*Cm2)*throttle - Cr2*v^2 - Cr1].*settings.delta_t


    return x_next
end
