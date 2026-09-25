# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Frame-aware quantity readers.
#
# Most orbital quantities are frame-dependent, and the older calc tags do not say
# which frame they are read in: they take whatever the spacecraft happens to
# carry. That is fine when everything is in one frame and silently wrong when
# it is not.
#
# These take any *subject*: a `Spacecraft`, an `AstroFrames.Coordinate`, or a
# user's own type implementing `state_of`/`frame_of`/`epoch_of`. An astronomer
# with a state vector and a script driving a vehicle write the same call.
#
# These are the lowercase readers, and the one home for the math. The older
# PascalCase tags such as `RAAN` do not call them: turning a zero-field tag into
# a struct that binds a subject and a frame would break every
# `OrbitCalc(sc, RAAN())` in the tree. The readers stand alone.
# =============================================================================

using AstroFrames: AbstractCoordinateSystem, Coordinate, frame_of
using AstroStates: CartesianState, KeplerianState, to_vector
using AstroUniverse: get_gravparam

"""
    position_vector(subject, cs::AbstractCoordinateSystem[, params])

Position of `subject` expressed in `cs`, in km.

# Arguments
- `subject`: a `Spacecraft`, an `AstroFrames.Coordinate`, or anything else
  implementing the subject interface. Its epoch and current frame come from it.
- `cs`: the coordinate system to read the position in.
- `params`: only for frames defined by data the frame subsystem cannot obtain
  on its own, such as an orbit-relative frame's `reference_state`.

# Returns
A three-element position vector in km.

# Notes
The frame is named rather than assumed. A position is meaningless without one,
and reading it out of whatever the subject happens to carry is how a quantity
ends up silently referred to the wrong frame.

# Examples
```julia
position_vector(sc, CoordinateSystem(earth, ITRF()))
position_vector(coord, CoordinateSystem(earth, MJ2000Ec()))
```
"""
function position_vector(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    return to_vector(CartesianState(Coordinate(subject, cs, params)))[1:3]
end

position_vector(subject, cs::AbstractCoordinateSystem) =
    position_vector(subject, cs, NamedTuple())

position_vector(subject) = position_vector(subject, frame_of(subject))

"""
    velocity_vector(subject, cs::AbstractCoordinateSystem[, params])

Velocity of `subject` expressed in `cs`, in km/s.

# Arguments
See [`position_vector`](@ref).

# Returns
A three-element velocity vector in km/s.

# Notes
Velocity is frame-dependent in a stronger sense than position: a rotating
target frame contributes a term of size `ω × r`, which at Earth's rate is
about 0.47 km/s at the equator. Reading a velocity in the wrong frame is an error of
that size, not a rounding difference.

# Examples
```julia
velocity_vector(sc, CoordinateSystem(earth, ITRF()))
```
"""
function velocity_vector(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    return to_vector(CartesianState(Coordinate(subject, cs, params)))[4:6]
end

velocity_vector(subject, cs::AbstractCoordinateSystem) =
    velocity_vector(subject, cs, NamedTuple())

velocity_vector(subject) = velocity_vector(subject, frame_of(subject))

"""
    raan(subject, cs::AbstractCoordinateSystem[, params])

Right ascension of the ascending node of `subject`, in `cs`, in radians.

# Arguments
See [`position_vector`](@ref).

# Returns
The RAAN in radians.

# Notes
RAAN is measured in the coordinate system's equatorial plane from its X axis,
so it is meaningless without naming the frame: the same orbit reads 112.31° in
`MJ2000Eq`, 270.10° in `MODEq`, 268.50° in `TODEq` and 82.46° in `ITRF`. This
is the clearest case for why a quantity reader takes a coordinate system: the
number changes and nothing about it says so.

The gravitational parameter comes from the coordinate system's origin, so it is
not a separate argument.

# Examples
```julia
raan(sc, CoordinateSystem(earth, MJ2000Eq()))
```
"""
function raan(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    μ = get_gravparam(cs.origin)
    return KeplerianState(CartesianState(Coordinate(subject, cs, params)), μ).raan
end

raan(subject, cs::AbstractCoordinateSystem) = raan(subject, cs, NamedTuple())

raan(subject) = raan(subject, frame_of(subject))
