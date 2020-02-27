#=
author: Linh Pham Thi
created on: 20.12.2019

This function uses the Ipopt and JuMP packages to solve a non-linear optimization problem. This problem is formulated in the classical MPC form. This is a state reference follower control.
=#
using Ipopt, JuMP, LinearAlgebra

function run_mpc(xr, x, ul, x_plan_prev, horizon_length)

    # time1=time()
    # Parameters
        delta_t = 0.1

        x_current = x

        Q = Diagonal(vec(ones(1,4)))         # Weighting matrix for states
        Q[3,3] = 0  # Set penalty of heading angle to zero due to now knowing reference
        R = Diagonal(vec(fill(4,2)))       # Weighting matrix for inputs

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

    # time2 = time()

    model = Model(with_optimizer(Ipopt.Optimizer, max_cpu_time=60.0))#, acceptable_tol= 10^-3)) #, print_level=0))

    # time3=time()

    # Variables
        @variable(model, x[1:4, 1:horizon_length])
        @variable(model, x_ref[1:4, 1:horizon_length])
        @variable(model, u[1:2, 1:horizon_length])
        @variable(model, u_last[1:2])

    # Warm start
        for i = 1:horizon_length
            for j=1:4
                set_start_value(x[j,i], x_plan_prev[j,2])
            end

            for j=1:2
                set_start_value(u[j,i], ul[j])
            end
        end

    # constraint the dummy variables to the input variables
    # we create these dummies to be used in the user-defined functions
    @constraint(model, x_ref .== xr)
    @constraint(model, u_last .== ul)

    # Initial conditions
        @constraint(model, [j=1:4], x[j,1] == x_current[j])

    # Bounds for input and speed
        @constraint(model, [i=1:horizon_length], 0<=(u[1,i])<=1)
        @constraint(model, [i=1:horizon_length], -1<=(u[2,i])<=1)

        @constraint(model, [i=1:horizon_length], x[4,i]>=0)

    # Dynamic constraints (bicycle model)

        @NLconstraint(model, [i=2:horizon_length],
                        x[1,i] == x[1, i-1] + (x[4, i-1] * cos(x[3, i-1]+C1 * u[2, i-1])*delta_t))
        @NLconstraint(model, [i=2:horizon_length],
                        x[2,i] == x[2, i-1] + (x[4, i-1] * sin(x[3, i-1]+C1 * u[2, i-1]))*delta_t)
        @NLconstraint(model, [i=2:horizon_length],
                        x[3,i] == x[3, i-1] + (C2 * u[2, i-1] * x[4, i-1])*delta_t)
        @NLconstraint(model, [i=2:horizon_length],
                x[4,i] == x[4, i-1] + (Cm1 * u[1,i-1] - Cm2 * u[1,i-1] * x[4,i-1] - Cr2 * x[4, i-1] * x[4, i-1])* delta_t - Cr1*delta_t)

    # TODO: add track boundary constaints

    # Setting objective function
    #
    # # METHOD 1: Use macro, can't feed gradients in, but it works
    #     @NLobjective(model, Min,
    #         sum(Q[j,j]*(x[j,i]-x_ref[j,i])^2 for i in 1:settings.horizon_length for j in 1:4) +
    #         sum(R[j,j]*(u[j,i]-u[j,i-1])^2 for i in 2:settings.horizon_length for j in 1:2) +
    #         sum(R[j,j]*(u[j,1]-u_last[j])^2 for j in 1:2))

    # METHOD 2: Register a function, use macro, doesn't work because the
    #           user-defined function must have scalar inputs

    # MG: HERE I AM ASSUMING DIAGONAL MATRICES R AND Q!!!!!
    function state_cost(x...)
        d = div(length(x), 2)
        return sum(Q[j, j] * (x[j] - x[j + d])^2 for j in 1:d)
    end

    function ∇state_cost(g, x...)
        d = div(length(x), 2)
        for j = 1:d
            g[j] = 2 * Q[j, j] * (x[j] - x[j + d])
        end
        for j = 1:d
            g[j + d] = -2 * Q[j, j] * (x[j] - x[j + d])
        end
    end

    function input_cost(u...)
        d = div(length(u), 2)
        return sum(R[j, j] * (u[j] - u[j + d])^2 for j in 1:d)
    end

    function ∇input_cost(g, u...)
        d = div(length(u), 2)
        for j = 1:d
            g[j] = 2 * R[j, j] * (u[j] - u[j + d])
        end
        for j = 1:d
            g[j + d] = -2 * R[j, j] * (u[j] - u[j + d])
        end
    end

    function input_last(u...)
        return R[1, 1] * (u[1] - u[3])^2 + R[2, 2] * (u[2] - u[4])^2
    end

    function ∇input_last(g, u...)
        g[1] = 2 * R[1, 1] * (u[1] - u[3])
        g[2] = 2 * R[2, 2] * (u[2] - u[4])
        g[3] = -2 * R[1, 1] * (u[1] - u[3])
        g[4] = -2 * R[2, 2] * (u[2] - u[4])
    end


    JuMP.register(model, :state_cost, 8, state_cost, ∇state_cost)
    JuMP.register(model, :input_cost, 4, input_cost, ∇input_cost)
    JuMP.register(model, :input_last, 4, input_last, ∇input_last)

    @NLobjective(model, Min,
    sum(state_cost(x[1, i], x[2, i], x[3, i], x[4, i], x_ref[1, i], x_ref[2, i], x_ref[3, i], x_ref[4, i]) for i in 1:settings.horizon_length) +
    sum(input_cost(u[1, i], u[2, i], u[1, i - 1], u[2, i - 1]) for i in 2:settings.horizon_length) +
    input_last(u[1, 1], u[2, 1], u_last[1], u_last[2]))


    # METHOD 3: Register a fuction, use the set_NL_objective function to set
    #           objective function instead, don't know the right syntax
        # JuMP.register(model, :my_objective_function, 22 + 10*horizon_length, state_ref_obj_func, autodiff=true)
        # JuMP.set_NL_objective(model, :Min, :my_objective_function,
        #                     [x[i,j] for i=1:4 for j=1:horizon_length])

    #

    # time4=time()

    JuMP.optimize!(model)

    # time5=time()

    # println("Time to set up parameters:", time2-time1, " s")
    # println("Time to set up model:", time3-time2, " s")
    # println("Time to set up model contraints and variables:", time4-time3, " s")
    # println("Time to solve:", time5-time4, " s")
    # println("")


    return JuMP.value.(x), JuMP.value.(u)[:,1]
end

# state_ref_obj_func(x, x_ref, u, u_last, Q, R, horizon_length) =
#     sum(a[5][j,j]*(a[1][j,i]-a[2][j,i])^2 for i in 1:horizon_length for j in 1:4) +
#     sum(a[6][j,j]*(a[3][j,i]-a[3][j,i-1])^2 for i in 2:horizon_length for j in 1:2) +
#     sum(a[6][j,j]*(a[3][j,1]-a[4][j])^2 for j in 1:2)
