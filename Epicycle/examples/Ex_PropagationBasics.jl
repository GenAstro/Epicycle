# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Propagation and Stopping Conditions
#'
#' Propagate a spacecraft to a duration, epoch, or orbital event, then draw the resulting
#' trajectory. A stopping condition names the quantity it watches, so the same quantities used to
#' inspect a state can also stop propagation.

using Epicycle
using EpicycleIO

#' ## Configuration
#'
#' Define the spacecraft, the forces acting on it, and the integrator the propagator uses. Drag and
#' solar radiation pressure each need geometry on the spacecraft as well as a force in the model:
#' the force says how the acceleration is computed, the geometry says what it acts on.

# Configure the spacecraft mass, drag area, and SRP area
sat = Spacecraft(state = CartesianState([5000.0, 5000.0, 0.0, -3.8, 3.8, 5.4]),
                 time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
                 coord_sys = CoordinateSystem(earth, ICRF()),
                 name = "Sat",
                 mass = 1000.0,
                 drag = SphericalDrag(c_d = 2.2, drag_area = 10.0),
                 srp = SphericalSRP(c_r = 1.8, srp_area = 10.0))

# Combine gravity with the Moon and Sun, atmospheric drag, and solar radiation
# pressure with a dual cone shadow
gravity = PointMassGravity(earth, (moon, sun))
drag = AtmosphericDrag(earth; model = Exponential())
srp = SolarRadiationPressure(earth; shadow = DualCone())
forces = ForceModel(gravity, drag, srp)

# Configure the propagator
integ = IntegratorConfig(Tsit5();
                         dt = 10.0,
                         reltol = 1e-9,
                         abstol = 1e-9)
prop = OrbitPropagator(forces, integ)

#' ## Propagate for a duration
#'
#' Duration is measured from the spacecraft's current epoch, in seconds or days.

# Propagate for one hour
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 3600.0))

# Propagate for half a day
propagate!(prop, sat, StopAt(sat, PropDurationDays(), 0.5))

#' ## Propagate to an epoch
#'
#' The epoch carries its own time scale, so a stop time in TDB stops a spacecraft kept in TAI.

# Stop at a fixed epoch
stop_time = Time("2015-09-22T12:00:00", TDB(), ISOT())
propagate!(prop, sat, StopAt(sat, stop_time))

#' ## Propagate to an orbital event
#'
#' A stopping condition is a quantity, a value, and the direction the quantity crosses it. Radius
#' dotted with velocity is zero at periapsis and apoapsis, and `direction` chooses between them:
#' the quantity increases through periapsis and decreases through apoapsis.

# Stop at periapsis and read the Keplerian state
propagate!(prop, sat, StopAt(position_dot_velocity, sat; equals = 0.0, direction = +1))
println(get_state(sat, Keplerian()))

# Stop at the ascending node, where z increases through zero
propagate!(prop, sat, StopAt(position_z, sat; equals = 0.0, direction = +1))
println(get_state(sat, Cartesian()))

# Propagate to apoapsis
propagate!(prop, sat, StopAt(position_dot_velocity, sat; equals = 0.0, direction = -1))
println(get_state(sat, Keplerian()))

# Stop at a radius of 7250 km and read the spherical RADEC state
propagate!(prop, sat, StopAt(position_magnitude, sat; equals = 7250.0))
println(get_state(sat, SphericalRADEC()))

#' ## Propagate backwards
#'
#' A negative duration or an earlier epoch runs the propagation backwards, with `direction = :infer`
#' taking the sign from the stopping condition.

# Propagate backwards for two hours
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), -7200.0); direction = :infer)

# Propagate backwards to an earlier epoch
stop_time_back = sat.time - 0.05
propagate!(prop, sat, StopAt(sat, stop_time_back); direction = :infer)

println(get_state(sat, Keplerian()))

#' ## Draw the trajectory
#'
#' `orbitview` reads the trajectory out of the spacecraft's own history and draws it on a globe,
#' which opens in a browser tab. The frame and the epochs come from the history, so the call takes
#' the spacecraft and a panel name.

# Draw the flown trajectory on the Cesium globe
orbitview("Trajectory", sat)
