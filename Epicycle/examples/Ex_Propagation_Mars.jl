# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Propagation About Mars
#'
#' Propagate a spacecraft in Mars orbit to periapsis. Both the spacecraft coordinate system and
#' the gravity model use Mars as their central body.

using Epicycle

#' ## Configuration
#'
#' Configure a spacecraft and point-mass propagator centered on Mars.

# Set the Mars-centered coordinate system
mars_icrf = CoordinateSystem(mars, ICRF())

# Initialize the spacecraft from Mars-centered Keplerian elements
elements = KeplerianState(6500.0,
                          0.3,
                          deg2rad(30),
                          deg2rad(145),
                          deg2rad(180),
                          deg2rad(0.01))

sat = Spacecraft(state = CartesianState(elements, mars.mu),
                 time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
                 coord_sys = mars_icrf,
                 name = "MarsSat")

# Configure unperturbed Mars gravity
gravity = PointMassGravity(mars, ())
forces = ForceModel(gravity)
integ = IntegratorConfig(Tsit5();
                         dt = 10.0,
                         reltol = 1e-9,
                         abstol = 1e-9)
prop = OrbitPropagator(forces, integ)

#' ## Propagate to periapsis

# Stop at periapsis and report the Mars-centered elements
propagate!(prop, sat, StopAt(position_dot_velocity, sat; equals = 0.0, direction = +1))
println(get_state(sat, Keplerian()))
