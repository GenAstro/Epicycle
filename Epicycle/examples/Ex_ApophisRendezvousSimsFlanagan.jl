# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Low-Thrust Earth-to-Apophis Rendezvous
#'
#' Design a 360-day rendezvous with asteroid 99942 Apophis using a one-newton engine.
#' The transfer uses Keplerian ephemerides for Earth and Apophis and minimizes
#' propellant with a `SimsFlanaganPhase`.
#'
#' Sims-Flanagan replaces continuous thrust with one impulse at the midpoint of
#' each segment. Half of the trajectory propagates forward from Earth and half
#' backward from Apophis; the two halves must meet at the match point. The solver
#' varies both control histories and the final mass while holding the initial
#' mass fixed.

using Epicycle
using LinearAlgebra
using Printf

#' ## Configuration
#'
#' Define the constants and initial conditions for the Sun, Earth, Apophis, and the spacecraft.

# Solar constants
const MU_SUN = 1.32712440018e11           # km^3/s^2
const AU = 1.495978707e8                  # km

# Earth osculating elements at departure
const EARTH_A = 1.00000011 * AU
const EARTH_E = 0.01671
const EARTH_I = 0.0
const EARTH_RAAN = 0.0
const EARTH_AOP = deg2rad(102.9372)
const EARTH_M = deg2rad(147.5124)

# Apophis osculating elements at departure
const APOPHIS_A = 0.9224 * AU
const APOPHIS_E = 0.1914
const APOPHIS_I = deg2rad(3.3393)
const APOPHIS_RAAN = deg2rad(203.9609)
const APOPHIS_AOP = deg2rad(126.7213)
const APOPHIS_M = deg2rad(178.8933)

# Spacecraft and propulsion
const M0 = 1500.0                         # kg
const ISP = 3000.0                        # s
const TMAX = 1.0e-3                       # kN, one newton
const G0 = 9.80665e-3                     # km/s^2

# Fixed transfer duration
const TOF = 360.0 * 86400.0               # s

#' ## Boundary Ephemerides
#'
#' Define the function to propagate Earth and Apophis orbits from their osculating elements.

# Advance the orbit from departure and return its Cartesian state and acceleration
function keplerian_ephemeris(a, e, i, raan, aop, m_departure)
    mean_motion = sqrt(MU_SUN / a^3)
    return function (t::Real)
        true_anomaly = mean_to_true_anomaly(m_departure + mean_motion * t, e)
        elements = KeplerianState(a, e, i, raan, aop, true_anomaly)
        state = CartesianState(elements, MU_SUN)
        r = collect(state.position)
        v = collect(state.velocity)
        return r, v, -(MU_SUN / norm(r)^3) .* r
    end
end

# Build the Earth and Apophis boundary ephemerides
earth_ephemeris = keplerian_ephemeris(EARTH_A, EARTH_E, EARTH_I,
                                      EARTH_RAAN, EARTH_AOP, EARTH_M)

apophis_ephemeris = keplerian_ephemeris(APOPHIS_A, APOPHIS_E, APOPHIS_I,
                                        APOPHIS_RAAN, APOPHIS_AOP, APOPHIS_M)

#' ## Configure the transfer
#'
#' The phase uses 40 segments, split evenly between the forward and backward
#' trajectories. Match-point scales put position, velocity, and mass residuals
#' on comparable numerical scales.
#'

const N_SEGMENTS = 40
const SMOOTHING = 0.015

_, v_departure, _ = earth_ephemeris(0.0)
_, v_arrival, _ = apophis_ephemeris(TOF)

phase = SimsFlanaganPhase(name = :earth_to_apophis,
                          transcription = SimsFlanagan(n_segments = N_SEGMENTS,
                                                       throttle_smoothing = SMOOTHING),
                          model = PropulsionModel(mu = MU_SUN,
                                                  Isp = ISP,
                                                  Tmax = TMAX,
                                                  g0 = G0),
                          ephemeris_left = earth_ephemeris,
                          ephemeris_right = apophis_ephemeris,
                          tspan = (0.0, TOF),
                          matchpoint_scale = [AU, AU, AU,
                                              norm(v_departure), norm(v_departure),
                                              norm(v_departure), M0])

# Initialize both control histories along their boundary velocities
u_departure = v_departure / norm(v_departure)
u_arrival = v_arrival / norm(v_arrival)

# Vary the two control histories and the boundary masses
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
     guess = 0.75 * M0,
     lower_bound = 100.0,
     upper_bound = M0,
     scale = M0,
     name = "mf")

# Limit the throttle magnitude to one in every segment
thrust_ball(c) = [dot(control(c), control(c))]
Constraint(thrust_ball, phase; lower_bound = 0.0, upper_bound = 1.0, at = Path())

#' ## The objective
#'
#' Propellant is proportional to throttle magnitude integrated over the transfer.
#' Include and epsilon (called SMOOTHING) to avoid singularities when the control
#' magnitude is near zero.

# Convert one segment at full throttle to propellant mass
const KG_PER_THROTTLE = TMAX * (TOF / N_SEGMENTS) / (G0 * ISP)

# Define the smoothed throttle magnitude 
throttle_magnitudes(u) = sqrt.(vec(sum(abs2, u, dims = 1)) .+ SMOOTHING^2)

# Define the objective and analytic partials for both control histories
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

#' ## Run the optimizer
#'
#' Solve the transfer problem using IPOPT

# Report the fixed mission inputs
@printf("m0    = %.1f kg\n", M0)
@printf("Tmax  = %.1f N\n", TMAX * 1000)
@printf("Isp   = %.0f s\n", ISP)
@printf("TOF   = %.1f days\n", TOF / 86400)

# Optimize the transfer
result = solve!(Sequence(phase); method = Optimize(max_iter = 500))

# Report delivered mass and the number of coast segments
mf = final_mass(phase)
magnitudes(u) = sqrt.(vec(sum(abs2, u, dims = 1)))
coasting = count(<(0.01), magnitudes(forward_control(phase))) +
           count(<(0.01), magnitudes(backward_control(phase)))

@printf("\nstatus         : %s\n", result.info)
@printf("delivered mass : %.3f kg\n", mf)
@printf("propellant     : %.3f kg\n", M0 - mf)
@printf("coasting       : %d of %d segments\n", coasting, N_SEGMENTS)
