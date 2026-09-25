# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Two Transcriptions in One Problem
#'
#' Fly the first half of a low-thrust transfer with Sims-Flanagan and the second half with
#' Hermite-Simpson collocation. A continuity constraint joins the phases into one optimization
#' problem that maximizes the final orbit radius.
#'
#' Neither phase knows what the other is. A transcription is a property of the arc it flies, so a
#' problem can use one where the control is a sequence of impulses and another where it is a
#' continuous history, and the node between them is an ordinary continuity constraint.

using Epicycle
using LinearAlgebra

#' ## Configuration
#'
#' Define the problem in canonical units: the gravitational parameter and the initial orbit radius
#' are one, and the spacecraft thrusts continuously for 3.32 time units.

# Set the model constants and phase durations
const MU = 1.0
const M0 = 1.0
const TMAX = 0.1405
const MDOT = 0.0749
const TF = 3.32
const TMID = TF / 2

# Set the departure state and state bounds
const Y0 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
const LO = [-5.0, -5.0, -5.0, -3.0, -3.0, -3.0]
const HI = [5.0, 5.0, 5.0, 3.0, 3.0, 3.0]

# Define the collocation state and control
struct RaiseState{T} <: AbstractState
    x::T
    y::T
    z::T
    vx::T
    vy::T
    vz::T
    m::T
end

struct RaiseControl{T} <: AbstractControl
    ux::T
    uy::T
    uz::T
end

# Collect the gravitational parameter, thrust, and mass flow
struct RaiseModel
    mu::Float64
    thrust::Float64
    mdot::Float64
end

model = RaiseModel(MU, TMAX, MDOT)

#' ## Equations of motion
#'
#' Two-body motion plus a thrust acceleration on the current mass. Both phases are built against
#' this one model.

# Evaluate the shared dynamics
function raise!(dy, y::RaiseState, u::RaiseControl, p, t, m::RaiseModel)
    r = sqrt(y.x^2 + y.y^2 + y.z^2)
    a = m.thrust / y.m
    dy[1] = y.vx
    dy[2] = y.vy
    dy[3] = y.vz
    dy[4] = -m.mu * y.x / r^3 + a * u.ux
    dy[5] = -m.mu * y.y / r^3 + a * u.uy
    dy[6] = -m.mu * y.z / r^3 + a * u.uz
    dy[7] = -m.mdot * sqrt(u.ux^2 + u.uy^2 + u.uz^2 + 1e-12)
end

#' ## The initial guess
#'
#' A guess has to satisfy the equations of motion, and Sims-Flanagan needs its two ends to be a
#' trajectory apart rather than merely plausible. Thrust along the velocity, integrated with RK4,
#' gives both.

# Fly the guess with thrust along the velocity
function flown(t_end; steps = 2000)
    f(y) = (r = norm(y[1:3]); v = norm(y[4:6]);
            vcat(y[4:6], -MU .* y[1:3] ./ r^3 .+ (TMAX / y[7]) .* y[4:6] ./ v, -MDOT))
    y, h = vcat(Y0, M0), t_end / steps
    for _ in 1:steps
        k1 = f(y)
        k2 = f(y .+ h / 2 .* k1)
        k3 = f(y .+ h / 2 .* k2)
        k4 = f(y .+ h .* k3)
        y = y .+ (h / 6) .* (k1 .+ 2k2 .+ 2k3 .+ k4)
    end
    return y
end

y_mid = flown(TMID)
y_end = flown(TF)

#' ## The Sims-Flanagan half
#'
#' Twenty segments, each with an impulse, propagated from both ends to a match point. Its variables
#' are the two end states, the masses and the per-segment controls. The control is a throttle, so
#' its magnitude stays inside the unit ball rather than on it, which has to be said outright
#' because the box on the components admits a magnitude of √3 at its corners.

# Build the Sims-Flanagan phase
first_half = SimsFlanaganPhase(name = :first_half,
                               transcription = SimsFlanagan(n_segments = 20),
                               model = PropulsionModel(mu = MU,
                                                       Isp = TMAX / MDOT,
                                                       Tmax = TMAX,
                                                       g0 = 1.0),
                               tspan = (0.0, TMID))

# Fix the departure state and vary the handoff state
Vary(state, first_half;
     guess = Y0,
     lower_bound = Y0,
     upper_bound = Y0)

Vary(final_state, first_half;
     guess = y_mid[1:6],
     lower_bound = LO,
     upper_bound = HI)

# Vary the segment controls and endpoint masses
Vary(forward_control, first_half;
     lower_bound = fill(-1.0, 3),
     upper_bound = fill(1.0, 3))

Vary(backward_control, first_half;
     lower_bound = fill(-1.0, 3),
     upper_bound = fill(1.0, 3))

Vary(initial_mass, first_half;
     lower_bound = M0,
     upper_bound = M0)

Vary(final_mass, first_half;
     guess = y_mid[7],
     lower_bound = 0.5,
     upper_bound = M0)

# Constrain the throttle magnitude to the unit ball, at every segment
throttle(c) = [dot(control(c), control(c))]
Constraint(throttle, first_half; lower_bound = 0.0, upper_bound = 1.0, at = Path())

#' ## The collocation half
#'
#' The same dynamics over the second half of the transfer. Its control is a thrust direction rather
#' than a throttle, so its magnitude is held at one.

# Build the Hermite-Simpson phase
second_half = CollocationPhase(name = :second_half,
                               transcription = HermiteSimpson(n_steps = 30),
                               dynamics = raise!,
                               model = model,
                               state = RaiseState,
                               control = RaiseControl,
                               tspan = (TMID, TF))

# Vary the state history and thrust direction
Vary(state, second_half;
     guess = hcat(y_mid, 0.5 .* (y_mid .+ y_end), y_end),
     lower_bound = vcat(LO, 0.5),
     upper_bound = vcat(HI, M0))

Vary(control, second_half;
     guess = repeat([0.0, 1.0, 0.0], 1, 3),
     lower_bound = fill(-1.0, 3),
     upper_bound = fill(1.0, 3))

# Constrain the thrust direction to unit magnitude, at every mesh point
unit_thrust(c) = [control(c).ux^2 + control(c).uy^2 + control(c).uz^2]
Constraint(unit_thrust, second_half; equals = 1.0, at = Path())

#' ## Join the phases
#'
#' `continuity` equates every state component and the time across the node between the phases,
#' which is what makes the second half start from the state the first produced.

# Join the phases with state, mass, and time continuity
seq = Sequence()
add_sequence!(seq, first_half)
add_sequence!(seq, second_half)
Constraint(continuity, Link(first_half, second_half))

# Arrive on a circular orbit, as high as the transfer can reach
final_radius(c) = (y = state(c); sqrt(y.x^2 + y.y^2 + y.z^2))
circular(c) = (y = state(c); [y.x * y.vx + y.y * y.vy + y.z * y.vz])

Constraint(circular, second_half; equals = 0.0, at = Final())
Objective(final_radius, second_half; sense = Max())

#' ## Solve the problem

# Solve and report the final orbit and handoff continuity
result = solve!(seq; method = Optimize(max_iter = 1000, tol = 1e-6, print_level = 5))

Y = state(second_half)
handoff = maximum(abs.(Y[:, 1] .- vcat(final_state(first_half), final_mass(first_half))))

println("status       : ", result.info)
println("final radius : ", round(norm(Y[1:3, end]), digits = 6))
println("final mass   : ", round(Y[7, end], digits = 6))
println("handoff jump : ", round(handoff, sigdigits = 3))
