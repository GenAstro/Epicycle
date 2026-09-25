# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Origin-coupling validation
#
# Some axes types are only physically meaningful at a specific origin (e.g.
# `ITRF` at Earth, `MoonPA` at the Moon). The `valid_origin` trait declares
# which combinations are allowed; `CoordinateSystem` construction checks it
# and errors with a domain-language message when the pair is invalid.
#
# The check is instance-based (not type-based) because all bodies in
# AstroUniverse share the same Julia type (`CelestialBody{Float64}`) and
# can only be distinguished at runtime — via NAIF ID for Earth/Moon axes
# and via the type parameter of `CelestialBodyFixed` for the planetary/Sun case.
#
# NAIF IDs used (from the standard SPICE catalogue):
#   Earth = 399, Moon = 301.
# =============================================================================

const _EARTH_NAIF_ID = 399
const _MOON_NAIF_ID  = 301

"""
    valid_origin(axes, origin) -> Bool

Whether `axes` makes sense about `origin`.

Some axes are meaningful only at a specific origin. For example, ITRF axes
require an Earth origin. Extension frames declare the same restriction by
adding a `valid_origin` method.

# Arguments
- `axes` — the axes instance.
- `origin` — the point the coordinate system would be centred on.

# Returns
`true` if the pair is a meaningful coordinate system, `false` otherwise.

# Notes
The default is `true`, so origin-independent frames require no additional
method.

Checked when the coordinate system is built, not when a transform runs, so the
mistake surfaces at the line that made it.

When a method returns `false`, `CoordinateSystem` raises an `ArgumentError` that
names the axes, the supplied origin, and the required origin.

# Example
```julia
AstroFrames.valid_origin(::PhobosFixed, origin) = origin.naifid == 401

CoordinateSystem(phobos, PhobosFixed())   # fine
CoordinateSystem(earth,  PhobosFixed())   # ArgumentError, naming the fix
```
"""
valid_origin(::AbstractAxes, ::Any) = true

# --- NAIF-ID helper ---------------------------------------------------------

@inline _naifid_of(o::AbstractPoint) =
    :naifid in propertynames(o) ? Int(getfield(o, :naifid)) : nothing

# --- Earth-restricted axes --------------------------------------------------

const _EARTH_RESTRICTED_AXES = (GCRF, CIRS, TIRS, ITRF, MODEq, TODEq, MODEc, TODEc, PEF, TEME)

for A in _EARTH_RESTRICTED_AXES
    @eval valid_origin(::$A, o) = (_naifid_of(o) == _EARTH_NAIF_ID)
end

# --- Moon-restricted axes ---------------------------------------------------

const _MOON_RESTRICTED_AXES = (MoonPA, MoonME)

for A in _MOON_RESTRICTED_AXES
    @eval valid_origin(::$A, o) = (_naifid_of(o) == _MOON_NAIF_ID)
end

# --- CelestialBodyFixed: origin's NAIF ID must match the type-parameter NAIF ID ---

valid_origin(::CelestialBodyFixed{N}, o) where {N} =
    (_naifid_of(o) == N)

# =============================================================================
# Error-message hints
# =============================================================================

_axes_family_hint(::Type{ITRF})  = "an Earth origin (use `earth`)"
_axes_family_hint(::Type{TIRS})  = "an Earth origin (use `earth`)"
_axes_family_hint(::Type{PEF})   = "an Earth origin (use `earth`)"
_axes_family_hint(::Type{TEME})  = "an Earth origin (use `earth`)"
_axes_family_hint(::Type{GCRF})  = "an Earth origin (use `earth`)"
_axes_family_hint(::Type{CIRS})  = "an Earth origin (use `earth`)"
_axes_family_hint(::Type{MODEq}) = "an Earth origin (use `earth`)"
_axes_family_hint(::Type{TODEq}) = "an Earth origin (use `earth`)"
_axes_family_hint(::Type{MODEc}) = "an Earth origin (use `earth`)"
_axes_family_hint(::Type{TODEc}) = "an Earth origin (use `earth`)"
_axes_family_hint(::Type{MoonPA}) = "a Moon origin (use `moon`)"
_axes_family_hint(::Type{MoonME}) = "a Moon origin (use `moon`)"
_axes_family_hint(::Type{<:CelestialBodyFixed{N}}) where {N} =
    "an origin whose NAIF ID matches the axes type parameter (NAIF ID = $(N))"
_axes_family_hint(::Type{<:AbstractAxes}) = "a specific origin type"
