# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # CCSDS OEM Ephemeris Files
#'
#' Read a CCSDS Orbit Ephemeris Message, write one, and confirm a written file holds the
#' trajectory it was written from. In this example we `read_ephemeris` an OEM that GMAT wrote,
#' propagate the same orbit from its first state and compare the two with `position_difference`.
#' Then we `write_ephemeris` the propagation, read the file back, and compare again.

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
#' GMAT wrote two hours of a low Earth orbit at 60 s steps, in EME2000 about the Earth. Each
#' segment of the file is a segment of the history, with its declared interpolation in the
#' segment's metadata.

# Read the file
gmat = read_ephemeris(joinpath(@__DIR__, "gmat_leo.oem"))
t_start, t_stop = span(gmat)
segment = gmat.segments[1]
println("epochs         : ", length(epochs(gmat)))
println("span           : ", t_start.isot, " to ", t_stop.isot, " UTC")
println("axes           : ", segment.coordinate_system.axes)
println("interpolation  : ", segment.metadata["interpolation"])

# Interpolate the state between two of the file's epochs
c = state_at(gmat, Time("2026-10-21T00:30:30", UTC(), ISOT()))
println("state at 00:30:30 : ", round.(to_vector(CartesianState(c)); digits = 6))

#' ## Propagate the same orbit
#'
#' GMAT's Earth gravitational parameter was set to Epicycle's, so the two propagations follow the
#' same orbit and differ only by the integrators.

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
println("GMAT OEM vs Epicycle, largest difference : ", round(maximum(d) * 1e3; sigdigits = 3), " m")

#' ## Write the propagation and read it back
#'
#' `step` writes a state every 60 s, interpolated from the propagation. The writer checks every
#' rule of the format before it opens the file, and a strict read accepts no correction, so a
#' file that reads back this way follows the standard exactly.

# Write the propagation as an OEM 3.0 file
path = joinpath(mktempdir(), "leo.oem")
write_ephemeris(path, sat;
                step = 60.0,
                object_id = "2026-001A",
                interpolation = Hermite(7))
foreach(println, first(readlines(path), 18))

# Read it back strictly and compare with the propagation
written = read_ephemeris(path; strict = true)
d = position_difference(written, sat)
println("epochs written : ", length(epochs(written)))
println("written OEM vs propagation, largest difference : ", round(maximum(d) * 1e9; sigdigits = 3), " µm")
println("interpolation  : ", written.segments[1].metadata["interpolation"])
