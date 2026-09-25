# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # The Spacecraft History
#'
#' Read a flown trajectory segment by segment. Every `propagate!` appends times, states, and the
#' recording coordinate system to the spacecraft history, so the mission script leaves its own
#' trajectory record.

using Epicycle
using LinearAlgebra

#' ## Configuration
#'
#' Configure a spacecraft, propagator, and transfer-orbit insertion maneuver.

# Initialize the spacecraft
sat = Spacecraft(state = CartesianState([3737.792, -4607.692, -2845.644, 5.411, 5.367, -1.566]),
                 time = Time("2000-01-01T11:59:28.000", UTC(), ISOT()),
                 name = "GeoSat-1")

# Configure the propagator
gravity = PointMassGravity(earth, ())
forces = ForceModel(gravity)
integ = IntegratorConfig(DP8();
                         abstol = 1e-12,
                         reltol = 1e-12,
                         dt = 60.0)
prop = OrbitPropagator(forces, integ)

# Define the transfer-orbit insertion maneuver
toi = ImpulsiveManeuver(axes = VNB(),
                        element1 = 2.518,
                        element2 = 0.0,
                        element3 = 0.0)

#' ## Fly the mission
#'
#' Coast, burn, coast. The two propagations leave one segment each; the maneuver between them is
#' what separates the segments.

# Propagate half a day, apply the maneuver, and propagate half a day more
propagate!(prop, sat, StopAt(sat, PropDurationDays(), 0.5))
maneuver!(sat, toi)
propagate!(prop, sat, StopAt(sat, PropDurationDays(), 0.5))

#' ## Read the history
#'
#' The history is indexed like a vector, and each segment carries its name.

# Inspect the recorded segments
println(sat.history)
println("Number of segments : ", length(sat.history))
println("First segment      : ", sat.history[1].name)
println("Last segment       : ", sat.history[end].name)

#' ## Read a segment
#'
#' A segment carries the epochs, the states, and the coordinate system they were recorded in.

# Inspect the first segment
seg = sat.history[1]
println("Name               : ", seg.name)
println("Number of points   : ", length(seg.times))
println("Origin             : ", seg.coordinate_system.origin.name)
println("Start              : ", seg.times[1])
println("End                : ", seg.times[end])
println("Duration (hours)   : ", (seg.times[end] - seg.times[1]) * 24)

#' ## Extract the trajectory
#'
#' Reading across the segments gives the whole flown trajectory, which is what a plot or a report
#' is built from.

# Collect positions and epochs across all segments
positions = [state.position for seg in sat.history for state in seg.states]
times = [time for seg in sat.history for time in seg.times]

# Split the positions into components and the epochs into modified Julian dates
x = [p[1] for p in positions]
y = [p[2] for p in positions]
z = [p[3] for p in positions]
times_mjd = [t.mjd for t in times]

println("Points in the whole trajectory: ", length(positions))
println("First position (km)           : ", positions[1])
println("Radius at the end (km)        : ", norm(positions[end]))
