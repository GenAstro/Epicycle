# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Earth to Mars, Low Thrust (Sims-Flanagan)
#'
#' Deliver as much mass as possible to Mars with a one-newton engine over a transfer time fixed at
#' the Hohmann value. A `SimsFlanaganPhase` connects circular, coplanar Earth and Mars ephemerides
#' while minimizing propellant use.
#'
#' Sims-Flanagan discretises an arc into segments, each carrying one impulse, propagated from both
#' ends to a match point in the middle. Its variables are not a state and a control history, which
#' is why they have their own names: `forward_control`, `backward_control`, `initial_mass` and
#' `final_mass`.
#'
#' The impulsive Hohmann transfer between the same orbits is the upper bound on delivered mass, and
#' the example prints both.

using Epicycle
using LinearAlgebra
using Printf

#' ## Configuration
#'
#' Use circular, coplanar orbits so the Earth and Mars ephemerides can be evaluated directly.

# Set the solar and planetary constants
const MU_SUN = 1.32712440018e11           # km^3/s^2
const AU = 1.495978707e8                  # km
const R_EARTH = 1.0 * AU
const R_MARS = 1.524 * AU

const VC_EARTH = sqrt(MU_SUN / R_EARTH)
const VC_MARS = sqrt(MU_SUN / R_MARS)
const OM_EARTH = VC_EARTH / R_EARTH
const OM_MARS = VC_MARS / R_MARS

# Set the spacecraft and engine properties
const M0 = 1500.0                         # kg
const ISP = 3000.0                        # s
const TMAX = 1.0e-3                       # kN, one newton
const G0 = 9.80665e-3                     # km/s^2

# Compute the Hohmann transfer time and the required initial Mars phase
const A_TRANSFER = (R_EARTH + R_MARS) / 2.0
const TOF = π * sqrt(A_TRANSFER^3 / MU_SUN)
const THETA_MARS0 = π - OM_MARS * TOF

#' ## The impulsive reference
#'
#' A two-burn Hohmann transfer between the same orbits delivers more mass than any low thrust
#' transfer can, because it spends its impulses where they are worth most. It is the upper bound
#' the answer is read against.

# Compute the impulsive upper bound
const DV1 = VC_EARTH * (sqrt(2 * R_MARS / (R_EARTH + R_MARS)) - 1)
const DV2 = VC_MARS * (1 - sqrt(2 * R_EARTH / (R_EARTH + R_MARS)))
const MF_IMPULSIVE = M0 * exp(-(DV1 + DV2) / (ISP * G0))

#' ## The ephemerides
#'
#' Each end of the phase follows a body. An ephemeris returns position, velocity and acceleration
#' at an epoch, which is all the transcription reads.

# Evaluate the Earth ephemeris
function earth_ephemeris(t::Real)
    θ = OM_EARTH * t
    r = R_EARTH * [cos(θ), sin(θ), 0.0]
    v = VC_EARTH * [-sin(θ), cos(θ), 0.0]
    a = -(MU_SUN / R_EARTH^2) * [cos(θ), sin(θ), 0.0]
    return r, v, a
end

# Evaluate Mars from the phase that places it at arrival
function mars_ephemeris(t::Real)
    θ = THETA_MARS0 + OM_MARS * t
    r = R_MARS * [cos(θ), sin(θ), 0.0]
    v = VC_MARS * [-sin(θ), cos(θ), 0.0]
    a = -(MU_SUN / R_MARS^2) * [cos(θ), sin(θ), 0.0]
    return r, v, a
end

#' ## Configure the transfer
#'
#' Sixty segments, thirty propagated forward from Earth and thirty backward from Mars. The match
#' point scale is what makes the two halves' disagreement comparable across position, velocity and
#' mass.

# Build the Sims-Flanagan phase
const N_SEGMENTS = 60
const SMOOTHING = 0.015

phase = SimsFlanaganPhase(name = :earth_to_mars,
                          transcription = SimsFlanagan(n_segments = N_SEGMENTS,
                                                       throttle_smoothing = SMOOTHING),
                          model = PropulsionModel(mu = MU_SUN,
                                                  Isp = ISP,
                                                  Tmax = TMAX,
                                                  g0 = G0),
                          ephemeris_left = earth_ephemeris,
                          ephemeris_right = mars_ephemeris,
                          tspan = (0.0, TOF),
                          matchpoint_scale = [R_EARTH, R_EARTH, R_EARTH,
                                              VC_EARTH, VC_EARTH, VC_EARTH, M0])

# Seed the throttles prograde at each end, which is the standard guess for an energy-raising transfer
_, v_departure, _ = earth_ephemeris(0.0)
_, v_arrival, _ = mars_ephemeris(TOF)
u_departure = v_departure / norm(v_departure)
u_arrival = v_arrival / norm(v_arrival)

# Vary both throttle blocks and the endpoint masses
Vary(forward_control, phase;
     guess = 0.5 * u_departure,
     lower_bound = fill(-2.0, 3),
     upper_bound = fill(2.0, 3),
     name = "u_fwd")

Vary(backward_control, phase;
     guess = 0.5 * u_arrival,
     lower_bound = fill(-2.0, 3),
     upper_bound = fill(2.0, 3),
     name = "u_bwd")

Vary(initial_mass, phase;
     lower_bound = M0,
     upper_bound = M0,
     name = "m0")

Vary(final_mass, phase;
     guess = 0.85 * M0,
     lower_bound = 100.0,
     upper_bound = M0,
     scale = M0,
     name = "mf")

# Constrain the throttle magnitude to the unit ball, at every segment
thrust_ball(c) = [dot(control(c), control(c))]
Constraint(thrust_ball, phase; lower_bound = 0.0, upper_bound = 1.0, at = Path())

#' ## The objective
#'
#' A segment burns propellant in proportion to its throttle magnitude, so the propellant is the
#' segment magnitudes summed and scaled to kilograms. The magnitude is smoothed by the same
#' `throttle_smoothing` the transcription burns with, so the objective prices a segment exactly as
#' the match point charges for it.

# Convert one unit of segment throttle to the propellant it burns
const KG_PER_THROTTLE = TMAX * (TOF / N_SEGMENTS) / (G0 * ISP)

# Compute the smoothed throttle magnitude of every segment in a block
throttle_magnitudes(u) = sqrt.(vec(sum(abs2, u, dims = 1)) .+ SMOOTHING^2)

# Define propellant use and its throttle partials
propellant(c) = KG_PER_THROTTLE * (sum(throttle_magnitudes(forward_control(c))) +
                                   sum(throttle_magnitudes(backward_control(c))))

@partial(propellant, forward_control) do c
    u = forward_control(c)
    KG_PER_THROTTLE .* vec(u ./ throttle_magnitudes(u)')
end

@partial(propellant, backward_control) do c
    u = backward_control(c)
    KG_PER_THROTTLE .* vec(u ./ throttle_magnitudes(u)')
end

Objective(propellant, phase; sense = Min())

#' ## Solve the transfer
#'
#' The problem solves at the solver's default settings, so the call names only the iteration limit.

# Report the transfer setup
@printf("m0    = %.1f kg\n", M0)
@printf("Tmax  = %.1f N\n", TMAX * 1000)
@printf("Isp   = %.0f s\n", ISP)
@printf("TOF   = %.2f days   (the Hohmann transfer time)\n", TOF / 86400)
@printf("impulsive Hohmann: dv1 = %.4f km/s, dv2 = %.4f km/s, mf = %.2f kg\n",
        DV1, DV2, MF_IMPULSIVE)

# Solve the transfer
result = solve!(Sequence(phase);
                method = Optimize(max_iter = 500, print_level = 5))

# Report the delivered mass against the impulsive bound
mf = final_mass(phase)
@printf("\nstatus         : %s\n", result.info)
@printf("delivered mass : %.3f kg   (impulsive bound %.2f kg)\n", mf, MF_IMPULSIVE)
@printf("propellant     : %.3f kg\n", M0 - mf)

# Report how much of the transfer coasts
magnitudes(u) = sqrt.(vec(sum(abs2, u, dims = 1)))
coasting = count(<(0.01), magnitudes(forward_control(phase))) +
           count(<(0.01), magnitudes(backward_control(phase)))
@printf("coasting       : %d of %d segments
", coasting, N_SEGMENTS)
