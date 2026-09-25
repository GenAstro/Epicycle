# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Targeting a Single Maneuver
#'
#' Size an impulsive maneuver so the following coast reaches a target apoapsis radius. The event
#' sequence varies the in-track burn component, propagates to apoapsis, and constrains its radius.

using Epicycle

#' ## Configuration
#'
#' Configure the spacecraft, propagator, and initial maneuver guess.

# Create the spacecraft
sat = Spacecraft()

# Configure the propagator
gravity = PointMassGravity(earth, (moon, sun))
forces = ForceModel(gravity)
integ = IntegratorConfig(Tsit5();
                         dt = 10.0,
                         reltol = 1e-9,
                         abstol = 1e-9)
prop = OrbitPropagator(forces, integ)

# Seed the VNB transfer-orbit insertion maneuver
toi = ImpulsiveManeuver(axes = VNB(),
                        element1 = 0.1)

#' ## Define the targeting problem
#'
#' Each `Event` is an action, its solver variables and its constraints. Only the in-track component
#' of the burn may move, which the bounds say by holding the other two at zero.

# Assemble the maneuver and coast events
seq = Sequence()
add_sequence!(seq,

    # Vary the in-track component of the burn
    Event(name = "TOI",
          event = () -> maneuver!(sat, toi),
          vars = [Vary(delta_v, toi;
                       lower_bound = [-10.0, 0.0, 0.0],
                       upper_bound = [10.0, 0.0, 0.0],
                       name = "toi")]),

    # Coast to apoapsis and constrain its radius
    Event(name = "Propagate to apoapsis",
          event = () -> propagate!(prop, sat,
                                   StopAt(position_dot_velocity, sat;
                                          equals = 0.0,
                                          direction = -1)),
          funcs = [Constraint(position_magnitude, sat; equals = 55000.0)]))

#' ## Solve the targeting problem

# Solve and report the sequence and solution
result = solve!(seq; method = Optimize(derivatives = :fd, print_level = 5))
report_sequence(seq)
report_solution(seq, result)
