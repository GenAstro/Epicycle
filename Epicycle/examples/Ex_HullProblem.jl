# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # The Hull Problem
#'
#' Minimize a Bolza objective with both terminal and running costs. A
#' Hermite-Simpson phase solves for a scalar state and its rate over a fixed
#' time interval.
#'
#' The terminal term penalizes the final state, while the path term integrates
#' control effort. The analytical solution ends at `x = 0.25` with a total cost
#' of `0.9375`.

using Epicycle

#' ## Problem Formulation
#'
#' The Bolza problem is
#'
#' ```math
#' \begin{aligned}
#' \underset{x(t),\,u(t)}{\operatorname{minimize}} \quad
#' & \frac{5}{2}x(1)^2 + \frac{1}{2}\int_0^1 u(t)^2\,dt \\
#' \text{subject to} \quad
#' & \dot{x} = u, \\
#' & x(0) = \frac{3}{2}.
#' \end{aligned}
#' ```
#'
#' where ``x`` is the scalar state and ``u`` is its rate. The analytical
#' solution is ``u(t)=-5/4`` and ``x(t)=3/2-5t/4``, which gives
#' ``x(1)=1/4`` and a total cost of ``15/16``.

#' ## Configuration
#'
#' Define a scalar state and a control equal to its time derivative. This simple
#' model keeps the example focused on the objective.

# Define the state and control variables
struct HullState{T} <: AbstractState
    x::T
end

struct HullControl{T} <: AbstractControl
    u::T
end

# Set the state rate equal to the control
hull_dynamics!(dy, y::HullState, u::HullControl, p, t, model) = (dy[1] = u.u)

# Supply the nonzero dynamics partial
@partial(hull_dynamics!, control) do dF, y, u, p, t, model
    dF[1, 1] = 1.0
end

#' ## Configure the Optimal Control Problem
#'
#' Fix the initial state and allow both the state and control histories to vary.
#' The optimizer balances final-state error against integrated control effort.

# Build the Hermite-Simpson phase
phase = CollocationPhase(name = :hull,
                         transcription = HermiteSimpson(n_steps = 20),
                         dynamics = hull_dynamics!,
                         state = HullState,
                         control = HullControl,
                         tspan = (0.0, 1.0))

# Vary the state and control histories
Vary(state, phase;
     guess = [1.5 0.25],
     lower_bound = [-10.0],
     upper_bound = [10.0])

Vary(control, phase;
     guess = reshape([-1.25 -1.25], 1, 2),
     lower_bound = [-10.0],
     upper_bound = [10.0])

# Define the initial-state constraint and both cost terms
initial_x(c) = [state(c).x]

terminal_cost(c) = 2.5 * state(c).x^2

running_cost(c) = 0.5 * control(c).u^2

# Supply analytic partials for the constraint and costs
@partial(initial_x, state) do c
    [1.0]
end

@partial(terminal_cost, state) do c
    [5.0 * state(c).x]
end

@partial(running_cost, control) do c
    [control(c).u]
end

# Fix the initial state
Constraint(initial_x, phase; equals = [1.5], at = Initial())

# Minimize terminal error and integrated control effort
Objective(terminal_cost, phase; sense = Min())
Objective(running_cost, phase; sense = Min(), at = Path())

#' ## Solve the Optimal Control Problem
#'
#' Check the supplied partials, solve the nonlinear program, and compare the
#' final state and objective with the analytical solution.

# Check the analytic partials against finite differences
check_partials(phase)

# Solve and report the final state and objective
result = solve!(Sequence(phase); method = Optimize(print_level = 5))
println("status    : ", result.info)
println("x_f       : ", round(get_final_state(phase).x, digits = 6), "   (analytical 0.25)")
println("objective : ", round(result.objective, digits = 6), "   (analytical 0.9375)")
