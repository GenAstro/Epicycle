# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Beta Angle, Eclipses and Sun–Moon Geometry
#'
#' Find the eclipses of a low Earth orbit and read its lighting geometry. In this example we create
#' a `Spacecraft`, read its beta angle, the angles among the Sun, the Earth, the Moon and the
#' spacecraft, and the fraction of the Sun hidden by the Earth and the Moon. Then we propagate it
#' to eclipse entry and exit with `StopAt` and plot the Earth's shadow over one day.

using Epicycle

#' ## Configuration

# Create a spacecraft on a 51.6 degree, 420 km circular orbit
orbit = KeplerianState(6798.137,
                       0.0005,
                       deg2rad(51.6),
                       deg2rad(30.0),
                       0.0,
                       0.0)
sat = Spacecraft(state = orbit,
                 time = Time("2026-06-21T00:00:00", UTC(), ISOT()),
                 name = "LEO")

# Configure a propagator with the Earth, Sun and Moon as point masses
gravity = PointMassGravity(earth, (sun, moon))
forces = ForceModel(gravity)
integ = IntegratorConfig(DP8();
                         dt = 30.0,
                         reltol = 1e-11,
                         abstol = 1e-11)
prop = OrbitPropagator(forces, integ)

#' ## Read the lighting geometry
#'
#' The beta angle is the elevation of the Sun above the orbit plane. It sets how long the eclipses
#' are: above about 70 degrees for this orbit there are none.

# Report the beta angles
println("beta angle               : ", round(rad2deg(beta_angle(sat)), digits = 4), " deg")
println("lunar beta angle         : ", round(rad2deg(lunar_beta_angle(sat)), digits = 4), " deg")
println("orbit normal to Sun      : ", round(rad2deg(orbit_normal_sun_angle(sat)), digits = 4), " deg")

# Report the angles among the Sun, the Earth, the Moon and the spacecraft
println("Sun-Earth-spacecraft     : ", round(rad2deg(sun_earth_spacecraft_angle(sat)), digits = 4), " deg")
println("Sun-spacecraft-Earth     : ", round(rad2deg(sun_spacecraft_earth_angle(sat)), digits = 4), " deg")
println("Sun-spacecraft-Moon      : ", round(rad2deg(sun_spacecraft_moon_angle(sat)), digits = 4), " deg")
println("Moon-Earth-Sun           : ", round(rad2deg(moon_earth_sun_angle(sat)), digits = 4), " deg")
println("Moon-Earth-spacecraft    : ", round(rad2deg(moon_earth_spacecraft_angle(sat)), digits = 4), " deg")
println("Moon-spacecraft-Earth    : ", round(rad2deg(moon_spacecraft_earth_angle(sat)), digits = 4), " deg")

# Report the shadow fractions: 0 in full sun, 1 in umbra
println("Earth shadow fraction    : ", earth_shadow_fraction(sat))
println("Moon shadow fraction     : ", moon_shadow_fraction(sat))

#' ## Find an eclipse
#'
#' The shadow fraction rises from 0 to 1 through the penumbra on the way into the Earth's shadow
#' and falls on the way out. Stopping where it crosses one half gives the middle of each
#' penumbra.

# Propagate to eclipse entry
propagate!(prop, sat, StopAt(earth_shadow_fraction, sat; equals = 0.5, direction = +1))
entry = sat.time
println("eclipse entry  : ", entry.utc.isot, " UTC")

# Propagate to eclipse exit
propagate!(prop, sat, StopAt(earth_shadow_fraction, sat; equals = 0.5, direction = -1))
exit_time = sat.time
println("eclipse exit   : ", exit_time.utc.isot, " UTC")
println("eclipse length : ", round((exit_time - entry) * 1440, digits = 2), " min")

#' ## Plot the shadow over one day

# Propagate one day
propagate!(prop, sat, StopAt(sat, PropDurationDays(), 1.0))

# Read the shadow fraction and beta angle along the flown trajectory
t, shadow, β = history(Calc(epoch, sat),
                       Calc(earth_shadow_fraction, sat),
                       Calc(beta_angle, sat))
hours = [(ti - t[1]) * 24 for ti in t]

# Plot the Earth shadow fraction
xyplot("Earth Shadow", hours, shadow; name = "LEO")
panel!("Earth Shadow"; xaxis_title = "hours", yaxis_title = "fraction of the Sun hidden")

# Report the range of the beta angle over the day
println("beta angle over the day : ", round(rad2deg(minimum(β)), digits = 3), " to ",
        round(rad2deg(maximum(β)), digits = 3), " deg")
