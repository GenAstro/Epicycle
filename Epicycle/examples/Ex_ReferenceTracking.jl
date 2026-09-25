# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Reference Tracking with a Path Bound
#'
#' Follow a reference signal at minimum cost while remaining inside a time-varying bound. An `LGL`
#' collocation phase solves the state and control histories over one period of the signal.
#'
#' The optimum is `x(t) = (2/3) sin t` at a cost of `π/3`, so the example has a closed-form value
#' to check against. The cost integrand and the path bound both depend on time explicitly, which
#' an optimal control problem posed on a trajectory rarely does.
#'
#' The `LGL` transcription is **Enterprise**, provided by the `EpicycleEnterprise` package.

using Epicycle
using EpicycleEnterprise

#' ## Configuration
#'
#' Define the state and the control. The state is the scalar being steered; the control is its rate.

# Set the final time to one reference period
const TF = 2π

# Define the state and control representations
struct TrackState{T} <: AbstractState
    x::T
end

struct TrackControl{T} <: AbstractControl
    u::T
end

#' ## Equations of motion
#'
#' The control is the rate of the state, so the dynamics carry no model of their own and the whole
#' problem is in the cost and the bound.

# Define the dynamics and control partial
track_dynamics!(dy, y::TrackState, u::TrackControl, p, t, model) = (dy[1] = u.u)

@partial(track_dynamics!, control) do dF, y, u, p, t, model
    dF[1, 1] = 1.0
end

#' ## The conditions and the cost
#'
#' The state returns to zero after one period, and along the way it stays inside an envelope whose
#' width follows the reference. The cost charges for straying from the reference and for control
#' effort, integrated over the phase.

# Define the endpoint state and its partial
position(c) = [state(c).x]

@partial(position, state) do c
    [1.0]
end

# Define the time-varying envelope and its state partial
envelope(c) = [(1.0 + 0.5 * sin(c.t)) * state(c).x^2 - 1.0]

@partial(envelope, state) do c
    [2.0 * (1.0 + 0.5 * sin(c.t)) * state(c).x]
end

# Define the cost integrand and its partials
tracking_error(c) = (state(c).x - sin(c.t))^2 + 0.5 * control(c).u^2

@partial(tracking_error, state) do c
    [2.0 * (state(c).x - sin(c.t))]
end

@partial(tracking_error, control) do c
    [control(c).u]
end

#' ## Configure the optimal control problem
#'
#' The state starts and ends at zero. The envelope is slack at the optimum, so it exercises a path
#' constraint without moving the answer.

# Build the LGL phase
phase = CollocationPhase(name = :tracking,
                         transcription = LGL(n_nodes = [24]),
                         dynamics = track_dynamics!,
                         state = TrackState,
                         control = TrackControl,
                         tspan = (0.0, TF))

# Vary the state and control histories
Vary(state, phase;
     guess = [0.0 0.0],
     lower_bound = [-5.0],
     upper_bound = [5.0])

Vary(control, phase;
     guess = reshape([0.0 0.0], 1, 2),
     lower_bound = [-5.0],
     upper_bound = [5.0])

# Pin both ends of the phase to zero
Constraint(position, phase; equals = [0.0], at = Initial())
Constraint(position, phase; equals = [0.0], at = Final())

# Hold the state inside the envelope at every node
Constraint(envelope, phase; upper_bound = [0.0], at = Path())

# Minimize the cost integrated over the phase
Objective(tracking_error, phase; sense = Min(), at = Path())

#' ## Solve the optimal control problem

# Check the declared partials against finite differences
check_partials(phase)

# Solve and compare the cost with the closed-form value
result = solve!(Sequence(phase); method = Optimize(print_level = 5))
println("status    : ", result.info)
println("cost      : ", round(result.objective, digits = 6), "   (analytic 1.047198)")
