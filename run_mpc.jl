#=
author: Linh Pham Thi
created on: 20.12.2019

This function uses the Ipopt and JuMP packages to solve a non-linear optimization problem. This problem is formulated in the classical MPC form. This is a state reference follower control.
=#
using Ipopt, JuMP, LinearAlgebra

function run_mpc(x_ref, x, u_last, horizon_length)

    # Parameters
        delta_t = 0.1

        x_current = x

        Q = Diagonal(vec(ones(1,4)))         # Weighting matrix for states
        Q[3,3] = 0  # Set penalty of heading angle to zero due to now knowing reference
        R = Diagonal(vec(fill(4,2)))         # Weighting matrix for inputs

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

    model = Model(with_optimizer(Ipopt.Optimizer, max_cpu_time=60.0))#, acceptable_tol= 10^-3)) #, print_level=0))

    # Variables
        @variable(model, x[1:4, 1:horizon_length])
        @variable(model, u[1:2, 1:horizon_length])


    # Warm start
        for i = 1:horizon_length
            for j=1:4
                set_start_value(x[j,i], x_ref[j])
            end

            for j=1:2
                set_start_value(u[j,i], u_last[j])
            end
        end


    # Initial conditions
        @constraint(model, [j=1:4], x[j,1] == x_current[j])

    # Bounds for input and speed
        @constraint(model, [i=1:horizon_length], 0<=(u[1,i])<=1)
        @constraint(model, [i=1:horizon_length], -1<=(u[2,i])<=1)

        @constraint(model, [i=1:horizon_length], x[4,i]>=0)
    # Dynamic constraints (bicycle model) +C1 * u[2, i-1]     x[4, i-1]

        @NLconstraint(model, [i=2:horizon_length],
                        x[1,i] == x[1, i-1] + (x[4, i-1] * cos(x[3, i-1]+C1 * u[2, i-1])*delta_t))
        @NLconstraint(model, [i=2:horizon_length],
                        x[2,i] == x[2, i-1] + (x[4, i-1] * sin(x[3, i-1]+C1 * u[2, i-1]))*delta_t)
        @NLconstraint(model, [i=2:horizon_length],
                        x[3,i] == x[3, i-1] + (C2 * u[2, i-1] * x[4, i-1])*delta_t)
#=
        @NLconstraint(model, [i=2:horizon_length],
                        x[4,i]*x[4,i-1] == x[4, i-1]*x[4,i-1] + (Cm1 * u[1,i-1] - Cm2 * u[1,i-1] * x[4,i-1] - Cr2 * x[4, i-1] * x[4, i-1])* x[4,i-1]*delta_t - Cr1*abs(x[4,i-1])*delta_t)
=#
        @NLconstraint(model, [i=2:horizon_length],
                x[4,i] == x[4, i-1] + (Cm1 * u[1,i-1] - Cm2 * u[1,i-1] * x[4,i-1] - Cr2 * x[4, i-1] * x[4, i-1])* delta_t - Cr1*delta_t)

    # TODO: add track boundary constaints

    @NLobjective(model, Min,
                sum(Q[j,j]*(x[j,i]-x_ref[j])^2 for i in 1:horizon_length for j in 1:4) +
                sum(R[j,j]*(u[j,i]-u[j,i-1])^2 for i in 2:horizon_length for j in 1:2) +
                sum(R[j,j]*(u[j,1]-u_last[j])^2 for j in 1:2))

    JuMP.optimize!(model)

    # println("MPC")
    # println("predicted x = ", JuMP.value.(x[:,2]))
    # println("u = ", JuMP.value.(u[:,1]))
    # println("objective function= ", JuMP.objective_value(model))
    # println(" ")

    return JuMP.value.(x), JuMP.value.(u)[:,1]
end
