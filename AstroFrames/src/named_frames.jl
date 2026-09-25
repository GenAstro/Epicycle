# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Frames under the names a GMAT user already knows.
#
# `inclination(sat, EarthMJ2000Ec)` rather than
# `inclination(sat, CoordinateSystem(earth, MJ2000Ec()))`. The pairing of an
# origin with axes is what a coordinate system is; these are the pairings
# common enough to deserve a name, spelled the way GMAT spells them so the
# mapping from an existing script is one for one.
#
# Anything else is still built the long way — these are a convenience, not a
# closed set.
# =============================================================================

"""
Earth-centred mean equator and equinox of J2000. GMAT's `EarthMJ2000Eq`.

# Example

```jldoctest
EarthMJ2000Eq

# output
CoordinateSystem:
  origin = Earth
  axes   = MJ2000Eq
```
"""
const EarthMJ2000Eq = CoordinateSystem(earth, MJ2000Eq())

"""
Earth-centred mean ecliptic and equinox of J2000. GMAT's `EarthMJ2000Ec`.

# Example

```jldoctest
EarthMJ2000Ec

# output
CoordinateSystem:
  origin = Earth
  axes   = MJ2000Ec
```
"""
const EarthMJ2000Ec = CoordinateSystem(earth, MJ2000Ec())

"""
Earth-centred Earth-fixed, ITRF. GMAT's `EarthFixed`.

# Example

```jldoctest
EarthFixed

# output
CoordinateSystem:
  origin = Earth
  axes   = ITRF
```
"""
const EarthFixed = CoordinateSystem(earth, ITRF())

"""
Earth-centred true of date equator. GMAT's `EarthTODEq`.

# Example

```jldoctest
EarthTODEq

# output
CoordinateSystem:
  origin = Earth
  axes   = TODEq
```
"""
const EarthTODEq = CoordinateSystem(earth, TODEq())

"""
Earth-centered ICRF, the inertial frame used by the ephemeris.

# Example

```jldoctest
EarthICRF

# output
CoordinateSystem:
  origin = Earth
  axes   = ICRF
```
"""
const EarthICRF = CoordinateSystem(earth, ICRF())

"""
Moon-centered mean Earth/mean rotation-axis coordinates for lunar mapping.

For gravity and dynamics use [`MoonPrincipalAxes`](@ref) instead; they differ
by about 875 m on the surface.

# Example

```jldoctest
MoonFixed

# output
CoordinateSystem:
  origin = Moon
  axes   = MoonME
```
"""
const MoonFixed = CoordinateSystem(moon, MoonME())

"""
Moon-centered principal-axis coordinates for lunar dynamics. See [`MoonFixed`](@ref).

# Example

```jldoctest
MoonPrincipalAxes

# output
CoordinateSystem:
  origin = Moon
  axes   = MoonPA
```
"""
const MoonPrincipalAxes = CoordinateSystem(moon, MoonPA())

"""
Sun-centred mean ecliptic and equinox of J2000.

# Example

```jldoctest
SunMJ2000Ec

# output
CoordinateSystem:
  origin = Sun
  axes   = MJ2000Ec
```
"""
const SunMJ2000Ec = CoordinateSystem(sun, MJ2000Ec())
