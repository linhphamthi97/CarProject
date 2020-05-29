#=
author: Linh Pham Thi
created on: 30.03.2020

This file contains two functions related to describing the track using the (s,d)
coordinate system. The first one is the B_spline function, which creates a B-spline
based on datapoints and a degree. This function uses the generate_N function, which
is a function to generate the B-spline coefficients.
=#
include("settings.jl")

using LinearAlgebra, CSV, DelimitedFiles

function B_spline(D, p, eval_points_no)
    #=
    This function creates a B-spline based on datapoints and a degree. It's made
    up of several stages.
    1. Select parameters t using the centripetal method
    2. Generate knot vector
    3. Generate matrix N
    4. Solve for matrix P
    5. Control points P and knot vector U and degree p determine an interpolating
       B-spline curve

    Inputs: D - datapoints
            p - degree of B-spline curve
            eval_points_no - number of evaluation points on the curve, usually set to 100
    Output: C - B-spline curve
    =#
    n = size(D,1) - 1   # n+1 is the number of data points

    # 1. Select parameters t
    global t = zeros(n + 1, 1)
    a = 0.5     # Select a= 1 for Chord Length Method, anything else for Centripetal Method
    L = sum(sqrt((D[i,j] - D[i-1,j])^2)^a for i in 2:size(D,1) for j in 1:2)

    t[1] = 0
    for k = 2:n
        global t[k] = sum(sqrt((D[i,j] - D[i-1,j])^2)^a for i in 2:k for j in 1:2)/L
    end
    t[n + 1] = 1

    # 2. Generate knot vector
    m = n + p + 1
    global U = zeros(m + 1, 1)

    for i in 1:m+1
        if i <= p+1
            U[i] = 0
        elseif i >= (m-p)+1
            U[i] = 1
        else
            # de Boor's averaging formula
            # U[i] = (1/p)* sum(t[k] for k in (i-p+1):(i))

            # Uniformly spaced method
            U[i] = (i-p-1)/(n-p+1)
        end
    end

    # 3. Generate matrix N
    N = zeros(n+1, n+1)
    for i in 1:size(t,1)
        N[i,:] = generate_N(n, p, m, t[i], U)
    end

    # 4. Solve for matrix P
    P = inv(N) * D

    # 5. Control points P and knot vector U and degree p determine an interpolating
    #       B-spline curve
    C = Array{Float64}(undef, settings.sample_number, 2)

    # For checking points (these points should coincide with the datapoints provided)
        # for i in 1:size(t,1)
        #     global C[i, :] = sum(N[j].*P[j,:] for j in 1:(n+1))
        # end

    i = 1
    for t in LinRange(0,1,eval_points_no)
        N = generate_N(n, p, m, t, U)
         global C[i,:] = sum(N[j].*P[j,:] for j in 1:(n+1))
         i += 1
     end

    return C, P, U
end

function generate_N(n, p, m, t, U)
    #=
    This function generates an array of basis functions for a B-spline of degree p
    evaluated at t.The detailed description behind this algorithm is on the website
    https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-curve-coef.html
    Essentially it is an implementation of the Cox-de Boor recursion formula.
    As mentioned on the website, this algorithm can be made more efficient (however
    I left it in its original form due to time constraints and for ease of understanding).

    Inputs: n - (no. of data points - 1)
            p - degree of B-spline curve
            m = n + p +1   - (m+1) is the number of knots
            t - a parameter
            U - knot vector
    Output: N - array of basis functions of dimensions (n+1)X1 evaluated at t
    =#

    # Initialise N to all zeros
    N = zeros(Float64, 1, n+1)

    # Rule out special cases
    if t == U[1]
        N[1] = 1.0
        return N
    elseif t == U[m+1]
        N[n+1] = 1.0
        return N
    end

    # Let t be in knot span u_k -> u_(k+1)
    if t == U[end]      # Special case
        global k = n
    else
        for a in 1:(size(U,1)-1)
            if U[a]<=t && U[a+1]>t
                global k = a - 1
            end
        end
    end

    # Degree 0 coefficient
    N[k+1] = 1.0

    # Rest of the coefficients
    for d in 1:p
        N[k-d+1] = ((U[k+2]-t)/(U[k+2]-U[k-d+2]))* N[k-d+2]

        for i in (k-d+1):(k-1)
            N[i+1] = ((t-U[i+1])/(U[i+d+1]-U[i+1])) * N[i+1] + ((U[i+d+2]-t)/(U[i+d+2]-U[i+2])) * N[i+2]
        end

        N[k+1] = ((t-U[k+1])/(U[k+d+1]-U[k+1])) * N[k+1]
    end

    return N
end

function length_calculator(P, n, p, U, t_begin, t_end)
    # Special Case
    if t_begin == t_end
        L = 0
        return L
    end

    m = n + p + 1

    if p == 1 # Rule out linear sections as that's easy to calculate based on end points
        N_begin = generate_N(n, p, m, t_begin, U)
        N_end = generate_N(n, p, m, t_end, U)
        C_begin = sum(N_begin[j].*P[j,:] for j in 1:(n+1))
        C_end = sum(N_end[j].*P[j,:] for j in 1:(n+1))

        L = sqrt((C_begin[1]-C_end[1])^2 + (C_begin[2]-C_end[2])^2)

    else
        t_vector = LinRange(t_begin, t_end, settings.integral_steps)
        C_deriv = zeros(settings.integral_steps, 2)

        for i in 1:settings.integral_steps
             t = t_vector[i]
        #     N = generate_N(n, (p - 1), m, t, U)
        #
        #     for k in 0:(n-1)
        #         C_deriv[i,:] += (N[k+1] * p / (U[k+p+2] - U[k+2])) *(P[k+2, :] - P[k+1, :])
        #     end

            C_deriv[i,:] = derivatives_calculator(P, n, p, U, t)[1:2]
            # println("C_deriv[i,:] for i = ", i, " is ", C_deriv[i,:])
        end

        L = sum(sqrt((C_deriv[i,1])^2 + (C_deriv[i,2])^2)*(t_end-t_begin)/settings.integral_steps for i in 1:settings.integral_steps)
    end
    return L
end

function derivatives_calculator(P, n, p, U, t)
    #=
    Inputs: P - control points
            n - number of datapoints - 1
            p - degree of B-spline curve
            U - knot vector
            t - parameter
    Output: C - Array of 1*4 containing the first (1:2) and second (3:4) derivatives
                of the B-spline curve at point parametrized by t
                    C[1]: dC/dx    C[2]: dC/dy      ... C[5]= curvature

    =#
    C_deriv = zeros(1,5)
    A = zeros(2,1)
    B = zeros(2,1)

    # Rule out linear sections as the derivatives are easy to calculate
    if p == 1
        N_begin = generate_N(n, p, n + p + 1, 0, U)
        N_end = generate_N(n, p, n + p + 1, 1, U)

        C_begin = sum(N_begin[j].*P[j,:] for j in 1:(n+1))
        C_end = sum(N_end[j].*P[j,:] for j in 1:(n+1))

        C_deriv[1] = (C_end[2]-C_begin[2])/(C_end[1]-C_begin[1])
        C_deriv[2] = 1/C_deriv[1]
        C_deriv[3:4] .= 0
    else
        # First derivative
        U1 = U[2:end-1]         # New knot sequence
        N = generate_N(n, (p - 1), size(U1, 1) -1, t, U1)

        for k in 0:(n-1)
            A += (P[k+2, :] - P[k+1, :]) .* (N[k+1]* p / (U[k+p+2] - U[k+2]))
        end
        C_deriv[1:2] = [A[2] A[1]]

        # Second derivative
        U2 = U[3:end-2]
        N = generate_N(n, (p - 2), size(U2, 2), t, U2)
        for k in 0:(n-2)
            B += (N[k+1] * (p-1) / (U[k+p+2] - U[k+3])) *
                            ((p*(P[k+3, :] - P[k+2, :])/(U[k+p+2] - U[k+2]))-(p*(P[k+2, :] - P[k+1, :])/(U[k+p+2] - U[k+2]))) # P'_i+1 - P'_i
        end
        C_deriv[3:4] = [B[2] B[1]]
    end

    # Curvature
    C_deriv[5] = abs(C_deriv[1]*C_deriv[4] - C_deriv[2]*C_deriv[3])/(C_deriv[1]^2 + C_deriv[2]^2)^(1.5)

    return C_deriv
end

# C = B_spline(D, 3, 100)
#
# plot(data_track[:,1], data_track[:,2], color=:black, xlim=(-2, 8), ylim=(-4, 2), label="Track data", legend=:bottomright)
# plot!(D[:,1], D[:,2], seriestype = :scatter, color=:blue, label="Data points", marker=(:hex, 5))
# plot!(C[:,1], C[:,2], color=:red, label="B-spline curve")
