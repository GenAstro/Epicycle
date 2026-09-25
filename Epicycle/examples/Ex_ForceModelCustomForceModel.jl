# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Low-Thrust Orbit Raising with a Custom Force Model
#'
#' Maximize the radius reached during a fixed-duration, low-thrust transfer. The
#' example defines a constant-thrust force, combines it with point-mass gravity,
#' and uses the resulting `ForceModel` as the dynamics of a Hermite-Simpson phase.
#'
#' The state includes spacecraft mass, and the control sets the thrust direction.
#' The departure state is fixed, the thrust direction remains normalized, and
#' the final orbit is constrained to be circular.

using Epicycle

# Import the acceleration function before extending it for the custom force
import AstroProp: accel_eval!

#' ## Configuration
#'
#' Define the phase state, thrust direction, and spacecraft. Mass is included in
#' the phase state because propellant use changes it along the trajectory.

# Define the phase state and control
struct RaiseState{T} <: AbstractState
    x::T
    y::T
    z::T
    vx::T
    vy::T
    vz::T
    m::T
end

struct ThrustDirection{T} <: AbstractControl
    ux::T
    uy::T
    uz::T
end

# Define the constant-thrust force
const G0 = 9.80665e-3

struct ConstantThrust <: OrbitODE
    thrust::Float64        # N
    isp::Float64           # s
end

# Add thrust acceleration and propellant mass flow to the state rates
function accel_eval!(model::ConstantThrust, t::Time, x̄, x̄̇, sc, params)
    m = x̄[7]
    u = params.control
    a = model.thrust * 1e-3 / m
    x̄̇[4] += a * u.ux
    x̄̇[5] += a * u.uy
    x̄̇[6] += a * u.uz
    throttle = sqrt(u.ux^2 + u.uy^2 + u.uz^2 + 1e-12)
    x̄̇[7] -= throttle * model.thrust * 1e-3 / (G0 * model.isp)
    return x̄̇
end

# Combine point-mass gravity with the custom thruster
thruster = ConstantThrust(200.0, 1500.0)
forces = ForceModel(PointMassGravity(earth, ()), thruster)

# Set the departure radius and fixed transfer duration
const R0 = 7000.0
const A_TARGET = 14000.0
const HALF_PERIOD = pi * sqrt(A_TARGET^3 / 398600.4415)

# Define the spacecraft used by the force model
sat = Spacecraft(state = CartesianState([R0, 0.0, 0.0, 0.0, 7.5460, 0.0]),
                 time = Time("2020-01-01T00:00:00.000", UTC(), ISOT()),
                 mass = 1000.0,
                 coord_sys = EarthICRF,
                 name = "raiser")

#' ## Build the Initial Guess
#'
#' Generate a dynamically consistent guess by integrating full thrust along the
#' velocity with RK4, then sampling the result at the collocation nodes. Starting
#' near a feasible trajectory leaves the optimizer to improve the transfer rather
#' than first repairing the dynamics.

# Propagate a full-thrust initial guess
function raise_guess(tf, n; steps = 4000)
    mu = 398600.4415
    thrust_km = thruster.thrust * 1e-3
    mdot = thrust_km / (G0 * thruster.isp)

    function f(y)
        r = y[1:3]
        v = y[4:6]
        rn = sqrt(r[1]^2 + r[2]^2 + r[3]^2)
        vn = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
        vcat(v, -mu / rn^3 .* r .+ (thrust_km / y[7]) .* (v ./ vn), -mdot)
    end

    y = [R0, 0.0, 0.0, 0.0, sqrt(mu / R0), 0.0, 1000.0]
    h = tf / steps
    Y = zeros(7, n)
    ts = range(0.0, tf; length = n)
    k, Y[:, 1] = 1, y
    for i in 1:steps
        k1 = f(y)
        k2 = f(y .+ h / 2 .* k1)
        k3 = f(y .+ h / 2 .* k2)
        k4 = f(y .+ h .* k3)
        y = y .+ (h / 6) .* (k1 .+ 2k2 .+ 2k3 .+ k4)
        t = i * h
        while k < n && t >= ts[k+1] - 1e-9
            k += 1
            Y[:, k] = y
        end
    end
    Y[:, n] = y
    return Y
end

#' ## Configure the Optimal Control Problem
#'
#' Use the combined force model as the phase dynamics and the spacecraft as the
#' object evaluated by those forces. Scale position, velocity, and mass by their
#' characteristic initial values.

# Set the position, velocity, and mass scales
const DU = R0
const VU = sqrt(398600.4415 / R0)
const MU_ = 1000.0

# Build the Hermite-Simpson phase
phase = CollocationPhase(name = :raise,
                         transcription = HermiteSimpson(n_steps = 40),
                         dynamics = forces,
                         model = sat,
                         state = RaiseState,
                         control = ThrustDirection,
                         tspan = (0.0, HALF_PERIOD))

# Vary the state history and thrust direction
Vary(state, phase;
     guess = raise_guess(HALF_PERIOD, 21),
     lower_bound = [-5e4, -5e4, -5e4, -20.0, -20.0, -20.0, 500.0],
     upper_bound = [5e4, 5e4, 5e4, 20.0, 20.0, 20.0, 1000.0],
     scale = [DU, DU, DU, VU, VU, VU, MU_])

Vary(control, phase;
     guess = repeat([0.0, 0.9, 0.0], 1, 21),
     lower_bound = [-1.5, -1.5, -1.5],
     upper_bound = [1.5, 1.5, 1.5])

# Define the quantities used by the constraints and objective
initial_state(c) = (y = state(c); [y.x, y.y, y.z, y.vx, y.vy, y.vz, y.m])

unit_direction(c) = [control(c).ux^2 + control(c).uy^2 + control(c).uz^2]

circular_arrival(c) = (y = state(c);
                       r = sqrt(y.x^2 + y.y^2 + y.z^2);
                       v2 = y.vx^2 + y.vy^2 + y.vz^2;
                       [y.x * y.vx + y.y * y.vy + y.z * y.vz, v2 - 398600.4415 / r])

final_radius_du(c) = (y = state(c); sqrt(y.x^2 + y.y^2 + y.z^2) / DU)

# Fix the departure, normalize thrust, and require a circular arrival
Constraint(initial_state, phase;
           equals = [R0, 0.0, 0.0, 0.0, 7.5460, 0.0, 1000.0],
           at = Initial(),
           scale = [DU, DU, DU, VU, VU, VU, MU_])

Constraint(unit_direction, phase; equals = 1.0, at = Path())

Constraint(circular_arrival, phase;
           equals = [0.0, 0.0],
           at = Final(),
           scale = [DU * VU, VU^2])

# Maximize final orbital radius
Objective(final_radius_du, phase; sense = Max())

#' ## Solve the Optimal Control Problem
#'
#' Solve the nonlinear program and report the final orbit characteristics

# Solve and report the final orbit
result = solve!(Sequence(phase); method = Optimize(tol = 1e-4, max_iter = 800, print_level = 5))
println("status       : ", result.info)
println("final radius : ", round(position_magnitude(subject_at(phase, Final())), digits = 3), " km")
println("final sma    : ", round(semi_major_axis(subject_at(phase, Final())), digits = 3), " km")
println("propellant   : ", round(1000.0 - get_final_state(phase).m, digits = 3), " kg")
println("time of flight: ", round(HALF_PERIOD, digits = 1), " s   (fixed)")
