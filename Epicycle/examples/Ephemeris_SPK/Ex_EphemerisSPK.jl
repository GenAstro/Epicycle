# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # JPL SPK Ephemeris Files
#'
#' Read a JPL SPK file, write one, and confirm a written file holds the trajectory it was written
#' from. In this example we `read_ephemeris` an SPK that GMAT wrote, propagate the same orbit from
#' its first state and compare the two with `position_difference`. Then we `write_ephemeris` the
#' propagation, read the file back, and compare again.

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
#' An SPK names its object by NAIF ID; this file holds one, so none needs to be given. Epochs are
#' SPICE ephemeris time, TDB seconds past J2000, and SPICE's J2000 frame reads as ICRF. A type 13
#' segment stores Hermite samples, which are read as stored.

# Read the file
gmat = read_ephemeris(joinpath(@__DIR__, "gmat_leo.bsp"))
t_start, t_stop = span(gmat)
segment = gmat.segments[1]
println("epochs         : ", length(epochs(gmat)))
println("span           : ", t_start.utc.isot, " to ", t_stop.utc.isot, " UTC")
println("object         : ", segment.metadata["object_id"])
println("axes           : ", segment.coordinate_system.axes)
println("interpolation  : ", segment.metadata["interpolation"])
foreach(d -> println("note           : ", d.message), diagnostics(gmat))

#' ## Propagate the same orbit

# Start a spacecraft from the file's first state
c0 = state_at(gmat, t_start)
sat = Spacecraft(state = CartesianState(c0),
                 time = t_start,
                 coord_sys = c0.coord_sys,
                 name = "LEO")

# Propagate for the two hours the file covers
propagate!(propagator, sat, StopAt(sat, PropDurationSeconds(), 7200.0))

# Compare the file with the propagation at the file's epochs
d = position_difference(gmat, sat)
println("GMAT SPK vs Epicycle, largest difference : ", round(maximum(d) * 1e3; sigdigits = 3), " m")

#' ## Write the propagation and read it back
#'
#' SPK needs the object's NAIF ID, negative for a spacecraft. States are written in ICRF axes,
#' the frame SPICE calls J2000, whatever the source's. A 64-bit ephemeris time resolves about
#' 0.1 µs this century, which bounds how closely the written file can match.

# Write the propagation as a Hermite (type 13) SPK, a state every 60 s
path = joinpath(mktempdir(), "leo.bsp")
write_ephemeris(path, sat;
                step = 60.0,
                naif_id = -900001,
                interpolation = Hermite(7))

# Read it back strictly and compare with the propagation
written = read_ephemeris(path; strict = true)
d = position_difference(written, sat)
println("epochs written : ", length(epochs(written)))
println("written SPK vs propagation, largest difference : ", round(maximum(d) * 1e6; sigdigits = 3), " mm")
println("interpolation  : ", written.segments[1].metadata["interpolation"])
