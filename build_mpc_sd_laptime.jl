#=
This function builds the MPC model through the JuMP interface. This file is for
the MPC model in the (s,d) coordinate system, with a state reference tracker.
=#
include("moving_ref_frame.jl")
include("track_s.jl")
using Ipopt, JuMP

struct MPCSolver

    model        #JuMP model
    xvar         #[4,N] array of JuMP variables for the state
    uvar         #[4,N] array of JuMP variables for the inputs
    xinit_param  #4 element array of JuMP parameters for the initial state x0
    uprev_param  #2 element array of JuMP parameters for the previous input u_{-1}

end

# Substitute curvature function
function gaussian(x, mu, sigma)
    return exp(-0.5*((x - mu)/sigma)^2)/(sigma*sqrt(2*pi))
end

function curvature_substitute(s)
    f = 2.4*gaussian(s,0,1) + gaussian(s,3.3,0.7) + 1.5*gaussian(s,6.8,0.7) +
        1.6*gaussian(s,9.8,0.3) + 1.5*gaussian(s,13.5,1) + 1.05*gaussian(s,15.7,0.3) +
        1.6*gaussian(s,18.9,0.3) + 1.6*gaussian(s,21.75,0.9) + 1.4*gaussian(s,25.5,1.5)
return f
end

function build_mpc(horizon_length)

    # Parameters
    Q_v = settings.Q_v
    R = settings.R_laptime       # Weighting matrix for inputs


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

    l_r = 0.05

    model = Model(with_optimizer(Ipopt.Optimizer,
    max_cpu_time=60.0,
    print_level=0,
    derivative_test = "none",
    hessian_approximation = "exact",
    print_options_documentation="no",
    print_user_options="no"))

    # Register curvature function
    register(model, :curvature_substitute, 1, curvature_substitute, autodiff = true)

    # Variables
    @variable(model, x[1:4, 1:horizon_length])
    @variable(model, u[1:2, 1:horizon_length])

    #Parameters
    @NLparameter(model, xinit[i=1:4] == 0.0)
    @NLparameter(model, uprev[i=1:2] == 0.0)

    #package model and its variables/parameters into a struct
    mpcsolver = MPCSolver(model,x,u,xinit,uprev)

    # Initial conditions
    @NLconstraint(model, [j=1:4], x[j,1] == xinit[j])

    # Bounds for inputs, speed, track boundary
    @constraint(model, [i=1:horizon_length], 0<=(u[1,i])<=1)
    @constraint(model, [i=1:horizon_length], -1<=(u[2,i])<=1)
    @constraint(model, [i=1:horizon_length], x[1,i]>=0)
    @constraint(model, [i=1:horizon_length], x[4,i]>=0)
    @constraint(model, [i=1:horizon_length], -settings.track_width/2<=x[2,i]<=settings.track_width/2)

    # Dynamic constraints (bicycle model)
    @NLconstraint(model, [i=2:horizon_length],
                    x[1,i] == x[1, i-1] + (x[4, i-1] * cos(x[3, i-1])/(1- curvature_substitute(x[1, i-1]) *x[2, i-1]))*settings.delta_t)
    @NLconstraint(model, [i=2:horizon_length], x[2,i] == x[2, i-1] + (x[4, i-1] * sin(x[3, i-1]))*settings.delta_t)
    @NLconstraint(model, [i=2:horizon_length], x[3,i] == x[3, i-1] + (x[4,i-1]*sin(atan(0.5*tan(C1*u[2, i-1])))/l_r -
                    (x[4, i-1] * cos(x[3, i-1]) * curvature_substitute(x[1, i-1]))/
                    (1- curvature_substitute(x[1, i-1]) *x[2, i-1]))*settings.delta_t)
    @NLconstraint(model, [i=2:horizon_length],
                    x[4,i] == x[4, i-1] + (Cm1 * u[1,i-1] - Cm2 * u[1,i-1] * x[4,i-1] - Cr2 * x[4, i-1] * x[4, i-1] - Cr1)*settings.delta_t)

    # Setting objective function
    @NLobjective(model, Min,
    -sum(x[1,i] for i in 1:settings.horizon_length) +
    sum(R[j,j]*(u[j,i]-u[j,i-1])^2 for i in 2:settings.horizon_length for j in 1:2) +
    sum(R[j,j]*(u[j,1]-uprev[j])^2 for j in 1:2)+
    Q_v*(x[4,end]))

    # JuMP.optimize!(model)

    return mpcsolver
end
