# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # GEO Transfer with an Event Sequence
#'
#' Target a geostationary transfer with three impulsive maneuvers and five
#' propagation events. Each maneuver event varies its delta-V, while propagation
#' events evaluate radius and inclination constraints at their stopping conditions.
#'
#' The bi-elliptic transfer first raises apoapsis above the mission orbit. A
#' mid-course correction lowers the inclination and places the next periapsis at
#' geostationary radius, where the final burn circularizes the orbit.

using Epicycle

#' ## Configuration
#'
#' Define the spacecraft, propagator, maneuver guesses, and stopping conditions.
#' Maneuver components are expressed in VNB axes.

# Set the initial spacecraft state
sat = Spacecraft(state = CartesianState([3737.792, -4607.692, -2845.644, 5.411, 5.367, -1.566]),
                 time = Time("2000-01-01T11:59:28.000", UTC(), ISOT()),
                 name = "GeoSat-1")

# Configure point-mass propagation
gravity = PointMassGravity(earth, ())
forces = ForceModel(gravity)
integ = IntegratorConfig(DP8();
                         abstol = 1e-12,
                         reltol = 1e-12,
                         dt = 60.0)
prop = OrbitPropagator(forces, integ)

# Set the initial maneuver guesses
toi = ImpulsiveManeuver(axes = VNB(),
                        element1 = 2.518,
                        element2 = 0.0,
                        element3 = 0.0)

mcc = ImpulsiveManeuver(axes = VNB(),
                        element1 = 0.559,
                        element2 = 0.588,
                        element3 = 0.0)

moi = ImpulsiveManeuver(axes = VNB(),
                        element1 = 0.282,
                        element2 = 0.0,
                        element3 = 0.0)

# Define the equator, apoapsis, and periapsis stopping conditions
z_crossing = StopAt(position_z, sat, EarthMJ2000Eq; equals = 0.0)
apoapsis = StopAt(position_dot_velocity, sat; equals = 0.0, direction = -1)
perigee = StopAt(position_dot_velocity, sat; equals = 0.0, direction = +1)

#' ## Define the Targeting Problem
#'
#' Events execute in the order shown below. Each `Vary` belongs to the maneuver
#' it changes, and each `Constraint` is evaluated after its event finishes.

# Assemble the event sequence
seq = Sequence()
add_sequence!(seq,

    # Coast to the first equatorial crossing
    Event(name = "Prop to Z 1",
          event = () -> propagate!(prop, sat, z_crossing)),

    # Apply transfer-orbit insertion and vary its tangential delta-V
    Event(name = "TOI",
          event = () -> maneuver!(sat, toi),
          vars = [Vary(delta_v, toi;
                       guess = [2.518, 0.0, 0.0],
                       lower_bound = [0.0, 0.0, 0.0],
                       upper_bound = [8.0, 0.0, 0.0],
                       name = "toi_v")]),

    # Coast to apoapsis and target an 85,000 km radius
    Event(name = "Prop to apoapsis",
          event = () -> propagate!(prop, sat, apoapsis),
          funcs = [Constraint(position_magnitude, sat; equals = 85000.0)]),

    # Coast to periapsis
    Event(name = "Prop to perigee 1",
          event = () -> propagate!(prop, sat, perigee)),

    # Continue to the next equatorial crossing
    Event(name = "Prop to Z 2",
          event = () -> propagate!(prop, sat, z_crossing)),

    # Apply the mid-course correction and vary its V and N components
    Event(name = "MCC",
          event = () -> maneuver!(sat, mcc),
          vars = [Vary(delta_v, mcc;
                       guess = [0.559, 0.588, 0.0],
                       lower_bound = [-1.0, -1.0, -0.001],
                       upper_bound = [4.0, 1.0, 0.001],
                       name = "mcc_vn")]),

    # Coast to periapsis and target its inclination and radius
    Event(name = "Prop to perigee 2",
          event = () -> propagate!(prop, sat, perigee),
          funcs = [Constraint(inclination, sat, EarthMJ2000Eq; equals = deg2rad(2.0)),
                   Constraint(position_magnitude, sat; equals = 42195.0)]),

    # Apply mission-orbit insertion and target the final semi-major axis
    Event(name = "MOI",
          event = () -> maneuver!(sat, moi),
          vars = [Vary(delta_v, moi;
                       guess = [0.282, 0.0, 0.0],
                       lower_bound = [-1.0, -0.001, -0.001],
                       upper_bound = [4.0, 0.001, 0.001],
                       name = "moi_v")],
          funcs = [Constraint(semi_major_axis, sat; equals = 42166.90)]))

#' ## Solve the Targeting Problem
#'
#' Solve with IPOPT using finite-difference partials, then report the event
#' sequence and converged maneuver values.

# Solve and report the transfer
result = solve!(seq; method = Optimize(max_iter = 1000,
                                       tol = 1e-6,
                                       derivatives = :fd,
                                       print_level = 5))
report_sequence(seq)
report_solution(seq, result)
