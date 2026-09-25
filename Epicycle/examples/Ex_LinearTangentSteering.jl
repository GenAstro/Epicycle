# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Linear Tangent Steering
#'
#' Steer a constant-acceleration vehicle to a target state in minimum time using a
#' `HermiteSimpson` collocation phase.
#'
#' The optimal steering angle obeys the linear tangent law, tan u = tan u₀ - c t. Solving that law
#' directly gives 0.554571 s, providing an independent check on the collocation result. Bryson and
#' Ho report 0.554 s.

using Epicycle
using LinearAlgebra

#' ## Problem Formulation
#'
#' The minimum-time problem is
#'
#' ```math
#' \begin{aligned}
#' \underset{x_1(t),\,x_2(t),\,x_3(t),\,x_4(t),\,u(t),\,t_f}
#' {\operatorname{minimize}} \quad & t_f \\
#' \text{subject to} \quad
#' & \dot{x}_1=x_2,\qquad \dot{x}_2=A\cos u, \\
#' & \dot{x}_3=x_4,\qquad \dot{x}_4=A\sin u, \\
#' & x_1(0)=x_2(0)=x_3(0)=x_4(0)=0, \\
#' & x_2(t_f)=45,\quad x_3(t_f)=5,\quad x_4(t_f)=0.
#' \end{aligned}
#' ```
#'
#' where ``x_1`` and ``x_3`` are position components, ``x_2`` and ``x_4``
#' are their velocities, ``A=100`` is the acceleration magnitude, and ``u``
#' is the steering angle. The final value of ``x_1`` is free.
#'
#' The optimal control satisfies the linear tangent law
#'
#' ```math
#' \tan u(t) = \tan u_0 - c t,
#' ```
#'
#' where ``u_0`` and ``c`` are constants determined by the terminal
#' conditions.

#' ## Configuration
#'
#' Use position and velocity in two axes as the state and the steering angle as the control.

# Set the constant acceleration magnitude
const A = 100.0

# Define the state and control representations
struct LTSState{T} <: AbstractState
    x1::T
    x2::T
    x3::T
    x4::T
end

struct LTSControl{T} <: AbstractControl
    u::T
end

#' ## Equations of motion
#'
#' The steering angle resolves the constant acceleration magnitude along the two axes.

# Evaluate the state rates
function lts_dynamics!(dy, y::LTSState, u::LTSControl, p, t, model)
    dy[1] = y.x2
    dy[2] = A * cos(u.u)
    dy[3] = y.x4
    dy[4] = A * sin(u.u)
end

# Supply the state partials
@partial(lts_dynamics!, state) do dF, y, u, p, t, model
    dF .= [0.0 1.0 0.0 0.0
           0.0 0.0 0.0 0.0
           0.0 0.0 0.0 1.0
           0.0 0.0 0.0 0.0]
end

# Supply the control partials
@partial(lts_dynamics!, control) do dF, y, u, p, t, model
    dF .= [0.0
           -A * sin(u.u)
           0.0
           A * cos(u.u)]
end

#' ## The initial guess
#'
#' Generate a dynamically consistent guess by integrating a linear steering sweep with RK4. The
#' sweep follows the expected shape of the linear tangent law.

# Fly the initial guess
const TF_GUESS = 0.6

u_guess(t) = 1.0 - 2.0 * t / TF_GUESS

function flown_guess(tf; steps = 400, knots = 5)
    f(y, t) = [y[2], A * cos(u_guess(t)), y[4], A * sin(u_guess(t))]
    h = tf / steps
    y = zeros(4)
    t = 0.0
    times = collect(range(0, tf; length = knots))
    Y = zeros(4, knots)
    next = 1
    for s in 0:steps
        while next <= knots && t >= times[next] - 1e-12
            Y[:, next] = y
            next += 1
        end
        s == steps && break
        k1 = f(y, t)
        k2 = f(y .+ h / 2 .* k1, t + h / 2)
        k3 = f(y .+ h / 2 .* k2, t + h / 2)
        k4 = f(y .+ h .* k3, t + h)
        y = y .+ h / 6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
        t += h
    end
    Y[:, end] = y
    return Y, times
end

Y_guess, t_guess = flown_guess(TF_GUESS)

#' ## Configure the optimal control problem
#'
#' The vehicle starts at rest at the origin and must arrive with a given velocity and crossrange.
#' The downrange position is free at the end, so only three of the four states are constrained.

# Build the Hermite-Simpson phase
phase = CollocationPhase(name = :linear_tangent_steering,
                         transcription = HermiteSimpson(n_steps = 60),
                         dynamics = lts_dynamics!,
                         state = LTSState,
                         control = LTSControl,
                         tspan = (0.0, TF_GUESS))

# Vary the state history, steering angle, and final time
Vary(state, phase;
     guess = Y_guess,
     lower_bound = [-10.0, -10.0, -10.0, -10.0],
     upper_bound = [100.0, 100.0, 100.0, 100.0],
     scale = [12.0, 45.0, 5.0, 10.0])

Vary(control, phase;
     guess = reshape(u_guess.(t_guess), 1, :),
     lower_bound = [-2.0],
     upper_bound = [2.0])

# A final time of zero is degenerate, so the lower bound stays clear of it
Vary(final_time, phase;
     guess = TF_GUESS,
     lower_bound = 0.1,
     upper_bound = 3.0)

# Define the boundary quantities and objective
launch(c) = [state(c).x1, state(c).x2, state(c).x3, state(c).x4]

arrival(c) = [state(c).x2, state(c).x3, state(c).x4]

duration(c) = final_time(c)

# Supply their analytic partials
@partial(launch, state) do c
    Matrix{Float64}(I, 4, 4)
end

@partial(arrival, state) do c
    [0.0 1.0 0.0 0.0
     0.0 0.0 1.0 0.0
     0.0 0.0 0.0 1.0]
end

@partial(duration, final_time) do c
    [1.0]
end

# Fix the launch and arrival conditions
Constraint(launch, phase; equals = [0.0, 0.0, 0.0, 0.0], at = Initial())
Constraint(arrival, phase; equals = [45.0, 5.0, 0.0], at = Final())

# Minimize the transfer time
Objective(duration, phase; sense = Min())

#' ## Solve the optimal control problem

# Check the declared partials against finite differences
check_partials(phase)

# Solve and compare the transfer time with the linear tangent law
result = solve!(Sequence(phase); method = Optimize(print_level = 5))
println("status : ", result.info)
println("tf     : ", round(get_final_time(phase), digits = 6),
        " s   (linear tangent law 0.554571 s)")
