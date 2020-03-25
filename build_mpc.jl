#=
This function builds the MPC model through the JuMP interface.
=#

using Ipopt, JuMP, LinearAlgebra

struct MPCSolver

    model        #JuMP model
    xvar         #[4,N] array of JuMP variables for the state
    uvar         #[4,N] array of JuMP variables for the inputs
    xref_param   #[4,N] array of JuMP parameters for the state reference trajectory
    xinit_param  #4 element array of JuMP parameters for the initial state x0
    uprev_param  #2 element array of JuMP parameters for the previous input u_{-1}

end


function build_mpc(horizon_length)

    # Parameters
    delta_t = 0.1

    Q = Diagonal(vec(ones(1,4)))         # Weighting matrix for states
    Q[3,3] = 0  # Set penalty of heading angle to zero due to not knowing reference
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

    model = Model(with_optimizer(Ipopt.Optimizer,
    max_cpu_time=60.0,
    print_level=0,
    derivative_test = "none",
    hessian_approximation = "exact",
    print_options_documentation="no",
    print_user_options="no"))


    # Variables
    @variable(model, x[1:4, 1:horizon_length])
    @variable(model, u[1:2, 1:horizon_length])

    #Parameters
    @NLparameter(model, xref[i=1:4, j=1:horizon_length] == 0.0)
    @NLparameter(model, xinit[i=1:4] == 0.0)
    @NLparameter(model, uprev[i=1:2] == 0.0)

    #package model and its variables/parameters into a struct
    mpcsolver = MPCSolver(model,x,u,xref,xinit,uprev)

    # Initial conditions
    @NLconstraint(model, [j=1:4], x[j,1] == xinit[j])

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
    @NLobjective(model, Min,
    sum(Q[j,j]*(x[j,i]-xref[j,i])^2 for i in 1:settings.horizon_length for j in 1:4) +
    sum(R[j,j]*(u[j,i]-u[j,i-1])^2 for i in 2:settings.horizon_length for j in 1:2) +
    sum(R[j,j]*(u[j,1]-uprev[j])^2 for j in 1:2))

    JuMP.optimize!(model)

    return mpcsolver
end
