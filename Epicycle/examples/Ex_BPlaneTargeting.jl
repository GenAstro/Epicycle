# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # B-Plane Targeting of an Earth Flyby
#'
#' Target the aim point of a hyperbolic Earth flyby with a trajectory correction maneuver. In this
#' example we create a `Spacecraft` inbound on a hyperbola about 190,000 km out, read its B-plane,
#' and build a `Sequence` that varies a correction maneuver so the B-plane at periapsis meets a
#' chosen B·T and B·R, with the Moon and Sun perturbing the approach.

using Epicycle

#' ## Configuration

# Create a spacecraft inbound on a hyperbola, about 190,000 km out
approach = KeplerianState(-20000.0,
                          1.4,
                          deg2rad(30.0),
                          deg2rad(60.0),
                          deg2rad(20.0),
                          deg2rad(-130.0))
sat = Spacecraft(state = approach,
                 time = Time("2026-09-01T00:00:00", UTC(), ISOT()),
                 name = "Flyby")

# Configure the propagator with the Earth, Moon and Sun as point masses
gravity = PointMassGravity(earth, (moon, sun))
forces = ForceModel(gravity)
integ = IntegratorConfig(DP8();
                         dt = 60.0,
                         reltol = 1e-11,
                         abstol = 1e-11)
prop = OrbitPropagator(forces, integ)

# Seed the trajectory correction maneuver
tcm = ImpulsiveManeuver(axes = VNB())

#' ## Read the B-plane
#'
#' T̂ lies along the incoming asymptote crossed with the coordinate system's Z axis, so B·T, B·R
#' and the B-plane angle depend on the frame they are read in. The B-vector magnitude, the impact
#' parameter, does not.

# Report the B-plane in equatorial and ecliptic axes
println("B·T, equatorial     : ", round(b_dot_t(sat), digits = 3), " km")
println("B·R, equatorial     : ", round(b_dot_r(sat), digits = 3), " km")
println("B-plane angle, equatorial : ", round(rad2deg(b_plane_angle(sat)), digits = 4), " deg")
println("B·T, ecliptic       : ", round(b_dot_t(sat, EarthMJ2000Ec), digits = 3), " km")
println("B·R, ecliptic       : ", round(b_dot_r(sat, EarthMJ2000Ec), digits = 3), " km")
println("B-vector magnitude  : ", round(b_magnitude(sat), digits = 3), " km")
println("hyperbolic excess speed : ", round(hyperbolic_excess_speed(sat), digits = 4), " km/s")
println("periapsis radius    : ", round(periapsis_radius(sat), digits = 3), " km")

#' ## Define the targeting problem
#'
#' The correction maneuver is free in all three components. The constraints are read at periapsis,
#' after the Moon and Sun have bent the approach.

seq = Sequence()
add_sequence!(seq,

    # Vary the correction maneuver
    Event(name = "TCM",
          event = () -> maneuver!(sat, tcm),
          vars = [Vary(delta_v, tcm;
                       lower_bound = [-0.2, -0.2, -0.2],
                       upper_bound = [0.2, 0.2, 0.2],
                       name = "tcm")]),

    # Coast to periapsis and constrain the B-plane there
    Event(name = "Propagate to periapsis",
          event = () -> propagate!(prop, sat,
                                   StopAt(position_dot_velocity, sat;
                                          equals = 0.0,
                                          direction = +1)),
          funcs = [Constraint(b_dot_t, sat; equals = 18500.0),
                   Constraint(b_dot_r, sat; equals = 3000.0)]))

#' ## Solve the targeting problem

# Solve and report the sequence and solution
result = solve!(seq; method = Optimize(derivatives = :fd, monitor = true))
report_sequence(seq)
report_solution(seq, result)

# Report the flyby the correction achieves
println("TCM delta-V         : ", round(delta_v_magnitude(tcm) * 1000, digits = 3), " m/s")
println("B·T at periapsis    : ", round(b_dot_t(sat), digits = 3), " km")
println("B·R at periapsis    : ", round(b_dot_r(sat), digits = 3), " km")
println("periapsis altitude  : ", round(altitude(sat), digits = 3), " km")
