# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Orbit Parameters from a Spacecraft State
#'
#' Read the orbit parameters of a Molniya orbit and of a departure hyperbola, set two of them, and
#' stop a propagation on one. In this example we create a `Spacecraft`, read its anomalies,
#' energy, period, apsis radii and speeds and its angular momentum, then set its mean anomaly and
#' apoapsis radius and propagate to apoapsis with `StopAt`.

using Epicycle

#' ## Configuration

# Create a spacecraft on a Molniya orbit
molniya = KeplerianState(26600.0,
                         0.74,
                         deg2rad(63.4),
                         deg2rad(40.0),
                         deg2rad(270.0),
                         deg2rad(30.0))
sat = Spacecraft(state = molniya,
                 time = Time("2026-03-01T00:00:00", UTC(), ISOT()),
                 name = "Molniya")

#' ## Read the orbit parameters
#'
#' Each parameter is read in the spacecraft's own coordinate system unless another is named.
#' Angles are in radians.

# Report the anomalies
println("true anomaly       : ", round(rad2deg(true_anomaly(sat)), digits = 4), " deg")
println("eccentric anomaly  : ", round(rad2deg(eccentric_anomaly(sat)), digits = 4), " deg")
println("mean anomaly       : ", round(rad2deg(mean_anomaly(sat)), digits = 4), " deg")

# Report the energy, period and mean motion
println("specific energy    : ", round(specific_energy(sat), digits = 4), " km²/s²")
println("C3                 : ", round(c3(sat), digits = 4), " km²/s²")
println("orbital period     : ", round(orbital_period(sat) / 3600, digits = 4), " h")
println("mean motion        : ", round(mean_motion(sat), sigdigits = 6), " rad/s")

# Report the apsides
println("periapsis radius   : ", round(periapsis_radius(sat), digits = 3), " km")
println("apoapsis radius    : ", round(apoapsis_radius(sat), digits = 3), " km")
println("periapsis speed    : ", round(periapsis_speed(sat), digits = 4), " km/s")
println("apoapsis speed     : ", round(apoapsis_speed(sat), digits = 4), " km/s")
println("semi-latus rectum  : ", round(semilatus_rectum(sat), digits = 3), " km")

# Report the angles measured in the orbit plane
println("argument of latitude   : ", round(rad2deg(argument_of_latitude(sat)), digits = 4), " deg")
println("longitude of periapsis : ", round(rad2deg(longitude_of_periapsis(sat)), digits = 4), " deg")
println("true longitude         : ", round(rad2deg(true_longitude(sat)), digits = 4), " deg")
println("descending node RA     : ", round(rad2deg(descending_node_right_ascension(sat)), digits = 4), " deg")

#' ## Read the vectors
#'
#' A vector is expressed in the axes of the coordinate system it is read in.

# Report the angular momentum, orbit normal and eccentricity vector in two frames
println("angular momentum, MJ2000Eq : ", round.(angular_momentum(sat, EarthMJ2000Eq), digits = 3), " km²/s")
println("angular momentum, MJ2000Ec : ", round.(angular_momentum(sat, EarthMJ2000Ec), digits = 3), " km²/s")
println("orbit normal               : ", round.(orbit_normal(sat), digits = 6))
println("eccentricity vector        : ", round.(eccentricity_vector(sat), digits = 6))
println("orbital angular velocity   : ", round.(orbital_angular_velocity(sat), sigdigits = 6), " rad/s")

#' ## Set orbit parameters
#'
#' A setter holds the rest of the representation the parameter belongs to: setting the mean
#' anomaly holds the other Keplerian elements, and setting the apoapsis radius holds the periapsis
#' radius and the orbit's orientation.

# Move the spacecraft to a mean anomaly of 90 degrees
mean_anomaly!(sat; to = deg2rad(90.0))
println("mean anomaly after setting : ", round(rad2deg(mean_anomaly(sat)), digits = 6), " deg")
println("true anomaly after setting : ", round(rad2deg(true_anomaly(sat)), digits = 4), " deg")

# Lower the apoapsis to 40,000 km
apoapsis_radius!(sat; to = 40000.0)
println("apoapsis radius  : ", round(apoapsis_radius(sat), digits = 3), " km")
println("periapsis radius : ", round(periapsis_radius(sat), digits = 3), " km")

#' ## Propagate to apoapsis
#'
#' The mean anomaly is π at apoapsis, so it serves as a stopping condition there. A stop at an
#' angle near its 0 / 2π wrap is not supported yet; stop at periapsis with `position_dot_velocity`.

# Configure a two-body propagator
forces = ForceModel(PointMassGravity(earth, ()))
integ = IntegratorConfig(DP8();
                         dt = 60.0,
                         reltol = 1e-12,
                         abstol = 1e-12)
prop = OrbitPropagator(forces, integ)

# Propagate to apoapsis and report the radius and speed there
propagate!(prop, sat, StopAt(mean_anomaly, sat; equals = Float64(π)))
println("radius at apoapsis : ", round(position_magnitude(sat), digits = 3), " km")
println("speed at apoapsis  : ", round(velocity_magnitude(sat), digits = 4), " km/s")
println("eccentric anomaly  : ", round(rad2deg(eccentric_anomaly(sat)), digits = 4), " deg")

#' ## Read a departure hyperbola
#'
#' On a hyperbola the hyperbolic anomaly replaces the eccentric anomaly, and the parameters of a
#' closed orbit, such as the period, throw an error saying why they do not apply.

# Create a spacecraft departing the Earth on a hyperbola
departure = Spacecraft(state = CartesianState([6778.0, 0.0, 0.0, 0.0, 11.6, 0.8]),
                       time = Time("2026-03-01T00:00:00", UTC(), ISOT()),
                       name = "Departure")

# Report the departure parameters
println("eccentricity            : ", round(eccentricity(departure), digits = 6))
println("hyperbolic anomaly      : ", round(hyperbolic_anomaly(departure), digits = 6), " rad")
println("mean anomaly            : ", round(mean_anomaly(departure), digits = 6), " rad")
println("hyperbolic excess speed : ", round(hyperbolic_excess_speed(departure), digits = 4), " km/s")
println("C3                      : ", round(c3(departure), digits = 4), " km²/s²")
println("periapsis speed         : ", round(periapsis_speed(departure), digits = 4), " km/s")

# A closed-orbit parameter on a hyperbola throws an ArgumentError
try
    orbital_period(departure)
catch err
    println(sprint(showerror, err))
end
