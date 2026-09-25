# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Hohmann Transfer
#'
#' Target a two-burn transfer from low Earth orbit to a circular orbit with a
#' radius of 45,000 km. The first maneuver raises apoapsis, the spacecraft coasts
#' to that apoapsis, and the second maneuver circularizes the orbit.

using Epicycle

#' ## Configuration
#'
#' Define the initial spacecraft, gravitational model, propagator, and maneuver
#' guesses. Both maneuvers are expressed in VNB axes.

# Set the initial spacecraft state
sat = Spacecraft(state = KeplerianState(7000.0, 0.001, 0.0, 0.0, 7.5, 1.0),
                 time = Time("2020-09-21T12:23:12", TAI(), ISOT()),
                 name = "Sat")

# Configure Earth, Moon, and Sun point-mass gravity
gravity = PointMassGravity(earth, (moon, sun))
forces = ForceModel(gravity)
integ = IntegratorConfig(DP8();
                         abstol = 1e-11,
                         reltol = 1e-11,
                         dt = 300.0)
prop = OrbitPropagator(forces, integ)

# Set the maneuver guesses
toi = ImpulsiveManeuver(axes = VNB(),
                        element1 = 0.1,
                        element2 = 0.2,
                        element3 = 0.3)

moi = ImpulsiveManeuver(axes = VNB(),
                        element1 = 0.4,
                        element2 = 0.5,
                        element3 = 0.6)

#' ## Define the Targeting Problem
#'
#' Vary the in-track component of each burn. At the end of the second maneuver,
#' constrain both orbital radius and eccentricity to define the target orbit.

# Assemble the maneuver-coast-maneuver sequence
seq = Sequence()
add_sequence!(seq,

    # Raise apoapsis with the first maneuver
    Event(name = "TOI",
          event = () -> maneuver!(sat, toi),
          vars = [Vary(delta_v, toi;
                       lower_bound = [0.0, 0.0, 0.0],
                       upper_bound = [2.5, 0.0, 0.0],
                       name = "toi")]),

    # Coast to apoapsis
    Event(name = "Coast to apoapsis",
          event = () -> propagate!(prop, sat,
                                   StopAt(position_dot_velocity, sat;
                                          equals = 0.0,
                                          direction = -1))),

    # Circularize and apply the terminal orbit constraints
    Event(name = "MOI",
          event = () -> maneuver!(sat, moi),
          vars = [Vary(delta_v, moi;
                       lower_bound = [0.0, 0.0, 0.0],
                       upper_bound = [3.0, 0.0, 0.0],
                       name = "moi")],
          funcs = [Constraint(position_magnitude, sat; equals = 45000.0),
                   Constraint(eccentricity, sat; equals = 0.0)]))

#' ## Solve the Targeting Problem
#'
#' Solve with IPOPT using finite-difference partials, then report the event
#' sequence and converged maneuver values.

# Solve and report the transfer
result = solve!(seq; method = Optimize(derivatives = :fd, print_level = 5))
report_sequence(seq)
report_solution(seq, result)
