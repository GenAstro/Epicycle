# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # STK Ephemeris Files
#'
#' Read an STK ephemeris (`.e`) file, write one, and confirm a written file holds the trajectory it
#' was written from. In this example we `read_ephemeris` a `.e` file that GMAT wrote, propagate
#' the same orbit from its first state and compare the two with `position_difference`. Then we
#' `write_ephemeris` the propagation with a maneuver in it, read the file back, and compare again.

using Epicycle

#' ## Configuration

# Configure the orbit propagator
forces = ForceModel(PointMassGravity(earth, ()))
integ = IntegratorConfig(DP8();
                         dt = 60.0,
                         reltol = 1e-13,
                         abstol = 1e-13)
propagator = OrbitPropagator(forces, integ)

#' ## Read the GMAT ephemeris
#'
#' STK times are seconds from `ScenarioEpoch`, in UTC. The file gives no `InterpolationSamplesM1`,
#' so STK's default, Lagrange degree 5, applies.

# Read the file
gmat = read_ephemeris(joinpath(@__DIR__, "gmat_leo.e"))
t_start, t_stop = span(gmat)
segment = gmat.segments[1]
println("epochs         : ", length(epochs(gmat)))
println("span           : ", t_start.utc.isot, " to ", t_stop.utc.isot, " UTC")
println("axes           : ", segment.coordinate_system.axes)
println("interpolation  : ", segment.metadata["interpolation"])

#' ## Propagate the same orbit
#'
#' GMAT computes the STK times from Modified Julian Dates, which carry about a microsecond, so
#' this file agrees with the propagation to millimetres rather than the micrometres of an exact
#' time column.

# Start a spacecraft from the file's first state
c0 = state_at(gmat, t_start)
sat = Spacecraft(state = CartesianState(c0),
                 time = t_start,
                 coord_sys = c0.coord_sys,
                 name = "LEO")

# Propagate for one hour, apply a maneuver, and propagate for the second hour
propagate!(propagator, sat, StopAt(sat, PropDurationSeconds(), 3600.0))
maneuver = ImpulsiveManeuver(axes = VNB(),
                             element1 = 0.01,
                             element2 = 0.0,
                             element3 = 0.0)
maneuver!(sat, maneuver)
propagate!(propagator, sat, StopAt(sat, PropDurationSeconds(), 3600.0))

# Compare the file with the propagation over the hour before the maneuver, the file's first 61 epochs
d = position_difference(gmat, sat)
println("GMAT .e vs Epicycle before the maneuver, largest difference : ", round(maximum(d[1:61]) * 1e3; sigdigits = 3), " m")

#' ## Write the propagation and read it back
#'
#' Each propagation is a segment, so the file carries a `SegmentBoundaryTimes` entry at the
#' maneuver, with the epoch repeated: the state before and the state after.

# Write the propagation, a state every 60 s
path = joinpath(mktempdir(), "leo.e")
write_ephemeris(path, sat; step = 60.0)
foreach(println, first(readlines(path), 16))

# Read it back strictly and compare with the propagation
written = read_ephemeris(path; strict = true)
d = position_difference(written, sat)
println("segments       : ", length(written.segments))
println("written .e vs propagation, largest difference : ", round(maximum(d) * 1e9; sigdigits = 3), " µm")

# The maneuver is a discontinuity in velocity at the shared epoch
t_burn = first(written.segments[2].times)
before = to_vector(written.segments[1].states[end])
after = to_vector(written.segments[2].states[1])
println("Δv at the maneuver : ", round(sqrt(sum(abs2, after[4:6] - before[4:6])) * 1e3; digits = 3), " m/s")
