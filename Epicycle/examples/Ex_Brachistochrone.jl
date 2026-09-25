# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # The Brachistochrone
#'
#' Find the fastest frictionless path between two points under constant gravity.
#' The problem uses Hermite-Simpson collocation to solve for the path, speed,
#' steering angle, and transfer time. Analytic partials are supplied for the
#' dynamics, constraints, and objective.
#'

using Epicycle

#' ## Problem Formulation
#'
#' The minimum-time problem is
#'
#' ```math
#' \begin{aligned}
#' \underset{x(t),\,y(t),\,v(t),\,\theta(t),\,t_f}{\operatorname{minimize}}
#' \quad & t_f \\
#' \text{subject to} \quad
#' & \dot{x} = v\sin\theta, \\
#' & \dot{y} = v\cos\theta, \\
#' & \dot{v} = g\cos\theta, \\
#' & x(0) = 0,\quad y(0) = 0,\quad v(0) = 0, \\
#' & x(t_f) = 2,\quad y(t_f) = 2, \\
#' & 0 \le v(t) \le 10,\quad
#'   -\frac{\pi}{2} \le \theta(t) \le \frac{\pi}{2}.
#' \end{aligned}
#' ```
#'
#' where ``x`` is horizontal position, ``y`` is downward vertical position,
#' ``v`` is speed, and ``\theta`` is the path angle measured from the
#' vertical.
#'
#' Without the inactive speed bound, the solution is a cycloid. Its
#' parameterization is
#'
#' ```math
#' x(\phi) = a(\phi - \sin\phi), \qquad
#' y(\phi) = a(1 - \cos\phi), \qquad
#' t(\phi) = \sqrt{\frac{a}{g}}\,\phi.
#' ```
#'
#' Choosing ``a`` and the final value of ``\phi`` to pass through ``(2,2)``
#' gives the reference transfer time of approximately ``0.8245`` seconds.

#' ## Configuration
#'
#' Define gravity, the bead state, and the path angle. The state contains
#' horizontal position, vertical position, and speed. The control is measured
#' from the vertical.

# Set the gravitational acceleration
struct BrachModel
    g::Float64
end

# Define the model with the gravitational acceleration
model = BrachModel(9.80665)

# Define the state and control variables
struct BrachState{T} <: AbstractState
    x::T
    y::T
    v::T
end

struct BrachControl{T} <: AbstractControl
    θ::T
end

#' ## Equations of Motion
#'
#' Resolve the bead velocity into horizontal and vertical components. Gravity
#' increases the speed through its component along the path.

# Evaluate the state rates
function brach_dynamics!(dy, y::BrachState, u::BrachControl, p, t, model)
    dy[1] = y.v * sin(u.θ)
    dy[2] = y.v * cos(u.θ)
    dy[3] = model.g * cos(u.θ)
end

# Supply the dynamics partial with respect to state
@partial(brach_dynamics!, state) do dF, y, u, p, t, model
    dF[1, 3] = sin(u.θ)
    dF[2, 3] = cos(u.θ)
end

# Supply the dynamics partial with respect to control
@partial(brach_dynamics!, control) do dF, y, u, p, t, model
    dF[1, 1] = y.v * cos(u.θ)
    dF[2, 1] = -y.v * sin(u.θ)
    dF[3, 1] = -model.g * sin(u.θ)
end

#' ## Configure the Optimal Control Problem
#'
#' The bead starts from rest at the origin and must reach `(2, 2)`. A
#' 20-step Hermite-Simpson mesh enforces the dynamics between optimization
#' nodes. The speed limit remains inactive at the solution.

# Build the Hermite-Simpson phase
phase = CollocationPhase(name = :brachistochrone,
                         transcription = HermiteSimpson(n_steps = 20),
                         dynamics = brach_dynamics!,
                         model = model,
                         state = BrachState,
                         control = BrachControl,
                         tspan = (0.0, 1.0))

# Vary the state history, path angle, and final time
Vary(state, phase;
     guess = [0.0 0.5 1.0 1.5 2.0
              0.0 0.5 1.0 1.5 2.0
              0.0 1.0 1.5 2.0 2.5],
     lower_bound = [-Inf, -Inf, 0.0],
     upper_bound = [Inf, Inf, Inf])

Vary(control, phase;
     guess = reshape([0.3, 0.5, 0.7, 0.9, 1.0], 1, 5),
     lower_bound = [-π / 2],
     upper_bound = [π / 2])

Vary(final_time, phase;
     lower_bound = 0.1,
     upper_bound = 10.0)

# Define the quantities used by the constraints and objective
start_state(c) = [state(c).x, state(c).y, state(c).v]

final_position(c) = [state(c).x, state(c).y]

speed(c) = [state(c).v]

duration(c) = final_time(c)

# Supply analytic partials for each quantity
@partial(start_state, state) do c
    [1.0 0.0 0.0
     0.0 1.0 0.0
     0.0 0.0 1.0]
end

@partial(final_position, state) do c
    [1.0 0.0 0.0
     0.0 1.0 0.0]
end

@partial(speed, state) do c
    [0.0 0.0 1.0]
end

@partial(duration, final_time) do c
    [1.0]
end

# Fix the endpoints and apply the speed limit by defining constraints
Constraint(start_state, phase; equals = [0.0, 0.0, 0.0], at = Initial())
Constraint(final_position, phase; equals = [2.0, 2.0], at = Final())
Constraint(speed, phase; upper_bound = 10.0, at = Path())
Constraint(speed, phase; upper_bound = 10.0, at = Final())

# Minimize transfer time
Objective(duration, phase; sense = Min())

#' ## Solve the Optimal Control Problem
#'
#' Check the supplied partials, solve the nonlinear program, and compare the
#' transfer time with the analytic cycloid solution.

# Check the analytic partials against finite differences
check_partials(phase)

# Solve and report the transfer time
result = solve!(Sequence(phase); method = Optimize())
println("status : ", result.info)
println("tf     : ", round(get_final_time(phase), digits = 6), " s   (cycloid 0.8245 s)")
