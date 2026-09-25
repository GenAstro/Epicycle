# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Parameter Identification from a Trajectory
#'
#' Recover a constant in the equations of motion from the trajectory it produced. An `LGL`
#' collocation phase varies the trajectory and parameter while minimizing the mismatch with the
#' data.
#'
#' The parameter is an optimization variable rather than a fixed property of the model, and the
#' phase has no control, so nothing steers the trajectory and the fit has one number to find. The
#' truth is `p = 0.1`, which puts the final state at 0.904837 and drives the cost to zero.
#'
#' The `LGL` transcription is **Enterprise**, provided by the `EpicycleEnterprise` package.

using Epicycle
using EpicycleEnterprise

#' ## Configuration
#'
#' Define the state. There is no control type, because the parameter is what the solver varies.

# Define the state representation
struct ParamIDState{T} <: AbstractState
    x::T
end

#' ## Equations of motion
#'
#' The state decays at a rate set by the parameter, so the partial with respect to the parameter is
#' the state itself.

# Define the dynamics and their state and parameter partials
pid_dynamics!(dy, y::ParamIDState, u, p, t, model) = (dy[1] = -p[1] * y.x)

@partial(pid_dynamics!, state) do dF, y, u, p, t, model
    dF[1, 1] = -p[1]
end

@partial(pid_dynamics!, parameter) do dF, y, u, p, t, model
    dF[1, 1] = -y.x
end

#' ## The condition and the cost
#'
#' The trajectory starts at one, and the cost is how far it sits from the data at each node.
#'
#' The cost integrand is left to automatic differentiation. It depends on time as well as on the
#' state, and writing its partial by hand buys nothing here, so the example declares partials where
#' they matter and lets the rest fall back. `check_partials` reports which ones did.

# Define the initial condition and its partial
initial_x(c) = [state(c).x]

@partial(initial_x, state) do c
    [1.0]
end

# Measure the trajectory mismatch
fit_error(c) = (state(c).x - exp(-0.1 * c.t))^2

#' ## Configure the optimal control problem
#'
#' The parameter is bounded well away from the truth on both sides, so the solver has to find it
#' rather than start on it.

# Build the LGL phase
phase = CollocationPhase(name = :param_id,
                         transcription = LGL(n_nodes = 20),
                         dynamics = pid_dynamics!,
                         state = ParamIDState,
                         tspan = (0.0, 1.0))

# Vary the state history and parameter
Vary(state, phase;
     guess = [1.0 exp(-0.1)],
     lower_bound = [-10.0],
     upper_bound = [10.0])

Vary(parameter, phase;
     guess = [0.5],
     lower_bound = [0.0],
     upper_bound = [2.0])

# Pin the trajectory to one at the initial time
Constraint(initial_x, phase; equals = [1.0], at = Initial())

# Minimize the mismatch integrated over the phase
Objective(fit_error, phase; sense = Min(), at = Path())

#' ## Solve the optimal control problem

# Check the declared partials against finite differences
check_partials(phase)

# Solve and compare the identified parameter with the truth
result = solve!(Sequence(phase); method = Optimize(print_level = 5))
println("status    : ", result.info)
println("parameter : ", round(get_param_value(phase)[1], digits = 6), "   (truth 0.1)")
println("x_f       : ", round(get_final_state(phase).x, digits = 6), "   (truth 0.904837)")
println("cost      : ", round(result.objective, digits = 8), "   (truth 0.0)")
