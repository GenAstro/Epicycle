# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # GEO Transfer with a Target Block
#'
#' Target the same three-burn geostationary transfer as the event-sequence
#' example, using a lightweight Domain-Specific Language (DSL) to create the DAG
#' using a`target!` block. Epicycle builds and solves the corresponding event sequence.
#' A `Vary` applies to the operation immediately after it. A `Constraint`
#' evaluates the state produced by the operation immediately before it.

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
integ = IntegratorConfig(Vern9();
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

#' ## Write the Flight Sequence
#'
#' Write the transfer in mission order and solve it with IPOPT using
#' finite-difference partials.

# Target the complete transfer
result = target!(method = Optimize(derivatives = :fd, print_level = 5)) do

    # Coast to the first equatorial crossing
    propagate!(prop, sat, z_crossing)

    # Apply transfer-orbit insertion and vary its tangential delta-V
    Vary(delta_v, toi; lower_bound = [0.0, 0.0, 0.0], upper_bound = [8.0, 0.0, 0.0])
    maneuver!(sat, toi)

    # Coast to apoapsis and target an 85,000 km radius
    propagate!(prop, sat, apoapsis)
    Constraint(position_magnitude, sat; equals = 85000.0)

    # Coast through periapsis to the next equatorial crossing
    propagate!(prop, sat, perigee)
    propagate!(prop, sat, z_crossing)

    # Apply the mid-course correction and vary its V and N components
    Vary(delta_v, mcc; lower_bound = [-1.0, -1.0, -0.001], upper_bound = [4.0, 1.0, 0.001])
    maneuver!(sat, mcc)

    # Coast to periapsis and target its inclination and radius
    propagate!(prop, sat, perigee)
    Constraint(inclination, sat, EarthMJ2000Eq; equals = deg2rad(2.0))
    Constraint(position_magnitude, sat; equals = 42195.0)

    # Apply mission-orbit insertion and target the final semi-major axis
    Vary(delta_v, moi; lower_bound = [-1.0, -0.001, -0.001], upper_bound = [4.0, 0.001, 0.001])
    maneuver!(sat, moi)
    Constraint(semi_major_axis, sat; equals = 42166.90)
end

#' ## Read the Solution
#'
#' The maneuver objects retain their solved delta-V values, and the spacecraft
#' retains the state produced by the final event.

# Report the solved burns and final orbit
println("status            : ", result.info)
println("TOI delta-v (km/s): ", round(delta_v(toi)[1], digits = 6))
println("MCC delta-v (km/s): ", round.(delta_v(mcc)[1:2], digits = 6))
println("MOI delta-v (km/s): ", round(delta_v(moi)[1], digits = 6))
println("final sma (km)    : ", round(semi_major_axis(sat), digits = 3))
