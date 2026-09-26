# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Getting Started
#'
#' Build a spacecraft, propagate it for a fixed duration, and read back the orbit it reached.
#' We describe the spacecraft with a `KeplerianState` and a `Time`, assemble a propagator from a
#' `ForceModel` and an `IntegratorConfig`, and end the propagation with a `StopAt` condition.

using Epicycle

#' ## The spacecraft
#'
#' A `Spacecraft` carries a state and the epoch that state belongs to. `KeplerianState` takes the
#' semi-major axis in kilometres followed by eccentricity, inclination, right ascension of the
#' ascending node, argument of periapsis and true anomaly, all in radians. The epoch names its own
#' time scale and format, here TAI written as an ISO 8601 string.

# Create a spacecraft from its orbit and epoch
sat = Spacecraft(
    state = KeplerianState(8000.0, 0.15, pi/4, pi/2, 0.0, pi/2),
    time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
    name = "sat",
)

#' ## The propagator
#'
#' A propagator is a force model and an integrator. `PointMassGravity` takes the central body
#' first and the bodies perturbing it after, so this is Earth gravity with the Moon and the Sun.
#' `IntegratorConfig` wraps any integrator from Julia's OrdinaryDiffEq: `dt` is the first step,
#' and the tolerances govern the adaptive steps after it.

# Configure Earth gravity with lunar and solar perturbations
gravity = PointMassGravity(earth, (moon, sun))
forces  = ForceModel(gravity)

# Configure the integrator and build the propagator
integ = IntegratorConfig(Tsit5(); dt = 10.0, reltol = 1e-9, abstol = 1e-9)
prop  = OrbitPropagator(forces, integ)

#' ## Propagate
#'
#' `StopAt` names the quantity that ends a propagation. `PropDurationSeconds` measures from the
#' spacecraft's current epoch, so this advances it 5000 seconds and leaves it there. The
#' spacecraft is updated in place, which is what the `!` marks, so reading its state afterwards
#' gives the orbit it arrived at. Angles go in as radians and are printed in degrees.

# Propagate for 5000 seconds
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 5000.0))

# Report the orbit it reached
println(get_state(sat, Keplerian()))
