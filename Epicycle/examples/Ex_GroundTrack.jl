# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Ground Track and Local Time of a Sun-Synchronous Orbit
#'
#' Follow a 700 km sun-synchronous orbit over the Earth and check that its local time of the
#' ascending node holds. In this example we create a `Spacecraft` at the ascending node, read its
#' altitude, geodetic latitude and longitude, sidereal and local time, and its J2 nodal precession
#' rate. Then we propagate it under J2, draw two orbits of ground track, and propagate thirty days
#' to compare the local time of the ascending node at the start and the end.

using Epicycle

#' ## Configuration

# Set the orbit: 700 km circular at the sun-synchronous inclination, starting at the ascending node
a = 7078.137
i = deg2rad(98.186)
v = sqrt(gravitational_parameter(earth) / a)
sat = Spacecraft(state = CartesianState([a, 0.0, 0.0, 0.0, v * cos(i), v * sin(i)]),
                 time = Time("2026-03-01T00:00:00", UTC(), ISOT()),
                 name = "SunSync")

# Configure the propagator with the Earth's J2 field
gravity = HarmonicGravity(earth;
                          degree = 2,
                          order = 0,
                          model = Zonal())
forces = ForceModel(gravity)
integ = IntegratorConfig(DP8();
                         dt = 60.0,
                         reltol = 1e-11,
                         abstol = 1e-11)
prop = OrbitPropagator(forces, integ)

#' ## Read the position over the Earth
#'
#' Altitude and latitude are geodetic, measured on the Earth's reference ellipsoid. Longitude is
#' east of Greenwich, and the velocity azimuth is the heading clockwise from north.

# Report the position relative to the ellipsoid
println("altitude           : ", round(altitude(sat), digits = 3), " km")
println("geodetic latitude  : ", round(rad2deg(geodetic_latitude(sat)), digits = 4), " deg")
println("longitude          : ", round(rad2deg(longitude(sat)), digits = 4), " deg")
println("velocity azimuth   : ", round(rad2deg(velocity_azimuth(sat)), digits = 4), " deg")
println("periapsis altitude : ", round(periapsis_altitude(sat), digits = 3), " km")
println("apoapsis altitude  : ", round(apoapsis_altitude(sat), digits = 3), " km")

# Report sidereal time and the mean geographic longitude
println("meridian hour angle       : ", round(rad2deg(meridian_hour_angle(sat)), digits = 4), " deg")
println("local sidereal time       : ", round(rad2deg(local_sidereal_time(sat)), digits = 4), " deg")
println("mean geographic longitude : ", round(rad2deg(mean_geographic_longitude(sat)), digits = 4), " deg")

#' ## Check the sun-synchronous condition
#'
#' A sun-synchronous orbit precesses eastward at the Sun's mean rate, one revolution per year.

# Compare the J2 nodal precession rate with the Sun's mean motion
sun_rate = 2π / (365.2422 * 86400)
println("nodal precession rate : ", round(nodal_precession_rate(sat), sigdigits = 6), " rad/s")
println("Sun's mean motion     : ", round(sun_rate, sigdigits = 6), " rad/s")

# Report the local time at the ascending node
println("mean local time : ", round(mean_local_time(sat), digits = 4), " h")
println("true local time : ", round(true_local_time(sat), digits = 4), " h")

#' ## Draw the ground track
#'
#' The history holds every propagated state, and `history` reads any quantity along it.

# Propagate two orbits
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 2 * orbital_period(sat)))

# Read the longitude, latitude and altitude along the flown trajectory
t, lon, lat, alt = history(Calc(epoch, sat),
                           Calc(longitude, sat),
                           Calc(geodetic_latitude, sat),
                           Calc(altitude, sat))

# Plot the ground track and the altitude
xyplot("Ground Track", rad2deg.(lon), rad2deg.(lat); mode = "markers", name = "SunSync")
panel!("Ground Track"; xaxis_title = "longitude (deg)", yaxis_title = "geodetic latitude (deg)")
minutes = [(ti - t[1]) * 1440 for ti in t]
xyplot("Altitude", minutes, alt; name = "SunSync")
panel!("Altitude"; xaxis_title = "minutes from start", yaxis_title = "altitude (km)")

#' ## Hold the local time of the ascending node
#'
#' The ascending node is where the geodetic latitude rises through zero. Thirty days later the
#' node has turned with the Sun, so its mean local time is nearly unchanged.

# Stop at the next ascending node and record its local time and RAAN
propagate!(prop, sat, StopAt(geodetic_latitude, sat; equals = 0.0, direction = +1))
mltan_start = mean_local_time(sat)
raan_start  = raan(sat)

# Propagate thirty days and stop at the ascending node again
propagate!(prop, sat, StopAt(sat, PropDurationDays(), 30.0))
propagate!(prop, sat, StopAt(geodetic_latitude, sat; equals = 0.0, direction = +1))
mltan_end = mean_local_time(sat)

# Report the node's drift against the Sun
println("mean local time of the node, start : ", round(mltan_start, digits = 4), " h")
println("mean local time of the node, end   : ", round(mltan_end, digits = 4), " h")
println("RAAN change over thirty days       : ", round(rad2deg(rem2pi(raan(sat) - raan_start, RoundNearest)), digits = 4), " deg")
println("Sun's mean motion over thirty days : ", round(rad2deg(sun_rate * 30 * 86400), digits = 4), " deg")
