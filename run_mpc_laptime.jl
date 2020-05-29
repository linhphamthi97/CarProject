#=
This function uses the Ipopt and JuMP packages to solve a non-linear optimization
 problem. This problem is formulated in the classical MPC form. This is the laptime
 minimizer.
=#

using Ipopt, JuMP, LinearAlgebra

function warm_start(mpcsolver, x_plan_prev, u_plan_prev, N)

    # Warm start state trajectory
    for i = 1:(N-1)
        for j=1:4
            set_start_value(mpcsolver.xvar[j,i], x_plan_prev[j,i+1])
        end
    end

    # Extend previous planning horizon by one
    xNp1 = simulation(x_plan_prev[:,N], u_plan_prev[:,N])
    for j=1:4
        set_start_value(mpcsolver.xvar[j,N], xNp1[j])
    end

    # warm start the input sequence
    for i = 1:N
        for j=1:2
            set_start_value(mpcsolver.uvar[j,i], u_plan_prev[j,1])
        end
    end

end

function run_mpc(mpcsolver, xinit, x_plan_prev, u_plan_prev, N)

    # warm start the x/u trajectories
    # warm_start(mpcsolver, x_plan_prev, u_plan_prev, N)

    #update parametric initial state
    for i = 1:4
        JuMP.set_value(mpcsolver.xinit_param[i], xinit[i])
    end

    #update parametric previous input
    for i = 1:2
        JuMP.set_value(mpcsolver.uprev_param[i], u_plan_prev[i,1])
    end

    JuMP.optimize!(mpcsolver.model)
    # println("Objective function value: ", objective_value(mpcsolver.model))
    # println("First obj func part: ", -sum(JuMP.value.(mpcsolver.xvar)[1,i] for i in 1:settings.horizon_length))
    # println("Third obj func part: ", sum(JuMP.value.(mpcsolver.xvar)[4,i] for i in 1:settings.horizon_length))
    return JuMP.value.(mpcsolver.xvar), JuMP.value.(mpcsolver.uvar)
end
