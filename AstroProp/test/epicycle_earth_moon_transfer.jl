# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The Epicycle side of gmat_earth_moon_transfer.script.

using Epicycle

earth_icrf = CoordinateSystem(earth, ICRF())
moon_icrf  = CoordinateSystem(moon, ICRF())

sat = Spacecraft(
    state     = CartesianState([2491.242350649993, -5494.004399970941, -2864.491425983717,
                                10.13590350345997, 3.454065939484507, 2.19038539314474]),
    time      = Time("2026-01-01T00:00:00", UTC(), ISOT()),
    coord_sys = earth_icrf,
)

integ      = IntegratorConfig(Vern9(); reltol = 1e-13, abstol = 1e-13)
earth_prop = OrbitPropagator(ForceModel(PointMassGravity(earth, (moon, sun))), integ)
moon_prop  = OrbitPropagator(ForceModel(PointMassGravity(moon, (earth, sun))), integ)

# Earth-centred to the Moon's sphere of influence
propagate!(earth_prop, sat, StopAt(position_magnitude, sat, moon_icrf; equals = 66100.0, direction = -1))

# Moon-centred to periapsis
propagate!(moon_prop, sat, StopAt(position_dot_velocity, sat, moon_icrf; equals = 0.0, direction = 1))

println(sat.time)
println(CartesianState(sat, moon_icrf))
