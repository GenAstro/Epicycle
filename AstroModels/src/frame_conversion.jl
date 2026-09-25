# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Expressing a spacecraft's state in another coordinate system.
#
# This is the call the frame subsystem exists to serve: *I have a spacecraft,
# give me its state somewhere else.*
#
# The conversion itself lives in AstroFrames, which owns it and does it for any
# subject. What lives here is the three methods that make a `Spacecraft` a
# subject, plus the `CartesianState` method. AstroFrames cannot name a
# `Spacecraft` without inverting the dependency direction, and AstroModels owns
# `Spacecraft`, so the methods live here and nothing is type piracy.
#
# It sits beside `get_state`, which is the same idea one axis over:
#
#     get_state(sc, Cartesian())     # change representation
#     CartesianState(sc, cs)         # change frame
# =============================================================================

using AstroFrames: Coordinate, CoordinateSystem, AbstractCoordinateSystem
import AstroFrames: state_of, frame_of, epoch_of
import AstroStates: CartesianState

# --- Spacecraft is a subject --------------------------------------------------
#
# Three lines, and every conversion AstroFrames offers works on a spacecraft,
# including any added later, without touching this file.

"""
    state_of(sc::Spacecraft)

The spacecraft's orbital state. See `AstroFrames.state_of`.
"""
state_of(sc::Spacecraft) = sc.state

"""
    frame_of(sc::Spacecraft)

The coordinate system the spacecraft's state is expressed in.
See `AstroFrames.frame_of`.
"""
frame_of(sc::Spacecraft) = sc.coord_sys

"""
    epoch_of(sc::Spacecraft)

The spacecraft's epoch. See `AstroFrames.epoch_of`.
"""
epoch_of(sc::Spacecraft) = sc.time

# --- The conversion -----------------------------------------------------------

"""
    CartesianState(sc::Spacecraft, cs::AbstractCoordinateSystem[, params]) -> CartesianState

The spacecraft's state expressed in `cs`.

Follows the conversion-constructor form used throughout `AstroStates`
(`CartesianState(kep, μ)`): the target type names the result, the source comes
first, and the extra data the conversion needs, here the coordinate system,
comes after.

The epoch is the spacecraft's own, so no epoch argument is needed.

# Arguments
- `sc`: the spacecraft, carrying its state, epoch, and current coordinate system.
- `cs`: the coordinate system to express the state in.
- `params`: evaluated numbers for an orbit-relative frame whose origin does not
  carry a reference orbit, such as `(; reference_state = x)` with `x` the reference
  position and velocity in ICRF, km and km/s. See `AstroFrames.axes_rotation`.

# Returns
A `CartesianState` in `cs`, with position in km and velocity in km/s.

# Notes
- Both halves of a frame change are applied: the axes rotation *and*, when the
  origin differs, the translation. Applying only the rotation is a silent error
  of the size of the origin separation.
- An origin change requires an ephemeris, and is therefore not differentiable.
  A pure rotation is.
- The same call works on an `AstroFrames.Coordinate`, for a script that has a
  state rather than a vehicle. Both go through one implementation.

# Examples
```julia
using AstroModels, AstroStates, AstroFrames, AstroUniverse
sc = Spacecraft()
CartesianState(sc, CoordinateSystem(earth, ITRF()))
CartesianState(sc, CoordinateSystem(moon,  MoonME()))

# Relative to another spacecraft: the chief is the origin and supplies the reference orbit
chief = Spacecraft(state = CartesianState([7000.0, -1.0, 0.0, 0.0, 7.5, 0.0]), name = "chief")
CartesianState(sc, CoordinateSystem(chief, RIC()))
```
"""
CartesianState(sc::Spacecraft, cs::AbstractCoordinateSystem, params::NamedTuple) =
    CartesianState(Coordinate(sc, cs, params))

CartesianState(sc::Spacecraft, cs::AbstractCoordinateSystem) =
    CartesianState(sc, cs, NamedTuple())
