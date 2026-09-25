# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

using Epicycle

# Create a spacecraft from its orbit and epoch
sat = Spacecraft(
    state=KeplerianState(8000.0,0.15,pi/4,pi/2,0.0,pi/2),
    time=Time("2015-09-21T12:23:12", TAI(), ISOT()),
    name = "sat",
)

# Configure Earth gravity with lunar and solar perturbations
gravity = PointMassGravity(earth,(moon,sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Propagate for 5000 seconds and report the final orbit
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 5000.0))
println(get_state(sat, Keplerian()))
