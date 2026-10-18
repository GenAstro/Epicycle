# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Nonsingular Keplerian and Delaunay Elements
#'
#' Describe a geostationary orbit with elements that stay defined where the Keplerian set does not.
#' In this example we create a `Spacecraft` on a near-circular, near-equatorial orbit, read its
#' Keplerian, nonsingular Keplerian and Delaunay elements, and propagate it to see which elements
#' change smoothly and which jump.

using Epicycle

#' ## Configuration

# Create a spacecraft on a near-geostationary orbit
geo = KeplerianState(42164.17,
                     2e-5,
                     deg2rad(0.01),
                     deg2rad(75.0),
                     deg2rad(20.0),
                     deg2rad(10.0))
sat = Spacecraft(state = geo,
                 time = Time("2026-09-01T00:00:00", UTC(), ISOT()),
                 name = "GEO")

# Configure the propagator with the Earth, Moon and Sun as point masses
gravity = PointMassGravity(earth, (moon, sun))
forces = ForceModel(gravity)
integ = IntegratorConfig(DP8();
                         dt = 300.0,
                         reltol = 1e-11,
                         abstol = 1e-11)
prop = OrbitPropagator(forces, integ)

#' ## Read the three element sets
#'
#' Near e = 0 and i = 0 the Keplerian RAAN, argument of periapsis and true anomaly are each poorly
#' defined, and only their sum carries the position. The nonsingular set holds the eccentricity and
#' inclination as vector components and the mean longitude as one angle.

# Report the Keplerian, nonsingular Keplerian and Delaunay elements
println(get_state(sat, Keplerian()))
println(get_state(sat, NonsingularKeplerian()))
println(get_state(sat, Delaunay()))

# Read single elements in a named coordinate system
println("mean longitude, MJ2000Eq : ", round(rad2deg(element(sat, NonsingularKeplerianState, EarthMJ2000Eq, :e5)), digits = 6), " deg")
println("Delaunay L               : ", round(element(sat, DelaunayState, EarthMJ2000Eq, :L), digits = 3), " km²/s")
println("Delaunay G               : ", round(element(sat, DelaunayState, EarthMJ2000Eq, :G), digits = 3), " km²/s")

#' ## Follow the elements over a week
#'
#' The Sun and Moon pull the orbit's eccentricity and inclination by small amounts. The Keplerian
#' node and periapsis swing through large angles as they do; the nonsingular elements move by the
#' small amounts the orbit actually changes.

# Propagate one week
propagate!(prop, sat, StopAt(sat, PropDurationDays(), 7.0))

# Read the Keplerian angles and two nonsingular elements along the flown trajectory
t, Ω, ω, e1, e3 = history(Calc(epoch, sat),
                          Calc(raan, sat),
                          Calc(argument_of_periapsis, sat),
                          Calc(element, sat, NonsingularKeplerianState, :e1),
                          Calc(element, sat, NonsingularKeplerianState, :e3))
days = [(ti - t[1]) for ti in t]

# Report the range each quantity covers in the week
println("RAAN range                 : ", round(rad2deg(maximum(Ω) - minimum(Ω)), digits = 3), " deg")
println("argument of periapsis range: ", round(rad2deg(maximum(ω) - minimum(ω)), digits = 3), " deg")
println("e1 range                   : ", round(maximum(e1) - minimum(e1), sigdigits = 3))
println("e3 range                   : ", round(maximum(e3) - minimum(e3), sigdigits = 3))

# Plot e1 and e3 over the week
xyplot("Nonsingular Elements", days, e1; name = "e1 = e cos(Ω + ω)")
xyplot!("Nonsingular Elements", days, e3; name = "e3 = sin(i/2) sin Ω")
panel!("Nonsingular Elements"; xaxis_title = "days", yaxis_title = "element value")
