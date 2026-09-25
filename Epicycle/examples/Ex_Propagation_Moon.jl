# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Propagation About the Moon
#'
#' Propagate a spacecraft in lunar orbit to periapsis, including perturbations from the Earth and
#' Sun. The central body changes in two places: the spacecraft coordinate system and the gravity
#' model.

using Epicycle

#' ## Configuration
#'
#' Configure a spacecraft and perturbed propagator centered on the Moon.

# Set the Moon-centered coordinate system
moon_icrf = CoordinateSystem(moon, ICRF())

# Initialize the spacecraft from Moon-centered Keplerian elements
elements = KeplerianState(5500.0,
                          0.38,
                          deg2rad(80),
                          deg2rad(145),
                          deg2rad(180),
                          deg2rad(0.01))

sat = Spacecraft(state = CartesianState(elements, moon.mu),
                 time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
                 coord_sys = moon_icrf,
                 name = "MoonSat",
                 cad_model = CADModel(file_path = joinpath(@__DIR__, "data", "DeepSpace1.obj"),
                                      scale = 100.0,
                                      visible = true))

# Configure lunar gravity with Earth and Sun perturbations
gravity = PointMassGravity(moon, (earth, sun))
forces = ForceModel(gravity)
integ = IntegratorConfig(Tsit5();
                         dt = 10.0,
                         reltol = 1e-9,
                         abstol = 1e-9)
prop = OrbitPropagator(forces, integ)

#' ## Propagate to periapsis
#'
#' Radius dotted with velocity increases through zero at periapsis, so the stopping condition is
#' the same one an Earth orbit uses.

# Stop at periapsis and report the Moon-centered elements
propagate!(prop, sat, StopAt(position_dot_velocity, sat; equals = 0.0, direction = +1))
println(get_state(sat, Keplerian()))
