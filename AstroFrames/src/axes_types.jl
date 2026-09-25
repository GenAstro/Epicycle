# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Axes types — community naming (ICRF/GCRF/ITRF style)
#
# Categories:
#
#   1. Body-fixed rotating axes            — restricted to specific origins
#      ITRF, TIRS, PEF                     (Earth only)
#      MoonPA, MoonME                      (Moon only; SPICE-backed)
#      CelestialBodyFixed{OT}                    (Sun, planets — origin encoded in type)
#
#   2. Origin-tied inertial axes           — restricted to Earth origin
#      GCRF, CIRS, MODEq, TODEq            (Earth-related model corrections)
#      MODEc, TODEc                        (Earth-orbit ecliptic-of-date)
#
#   3. Origin-agnostic inertial axes       — any origin allowed
#      ICRF, MJ2000Eq, MJ2000Ec
#
#   4. Orbit-relative axes                 — spacecraft-referenced
#      VNB                                 (existing)
#
# Origin-coupling rules are declared alongside each type as a
# `_valid_origin(::Type{A}, ::Type{O}) = true|false` method. Defaults to
# `true` (any origin allowed); overridden per-type as needed.
# =============================================================================

# --- Modern (IAU-2006 family) -----------------------------------------------

"""
    ICRF()

International Celestial Reference System axes.

Origin-agnostic inertial. Standards-strict usage pairs ICRF with the Solar
System Barycenter as origin; this framework permits any origin.

Shares its orientation exactly with `GCRF`: the GCRS is defined as
kinematically non-rotating with respect to the ICRS, so no rotation separates
them. They differ in origin, not in axes.

See also: [`GCRF`](@ref).

# Example

```jldoctest
using AstroUniverse: earth
CoordinateSystem(earth, ICRF())

# output
CoordinateSystem:
  origin = Earth
  axes   = ICRF
```
"""
struct ICRF <: AbstractAxes end

"""
    GCRF()

Geocentric Celestial Reference System axes.

Origin-tied inertial, and the geocentric frame the Earth chain is built on.
Requires an Earth origin.

Its orientation is identical to `ICRF`; the GCRS is kinematically
non-rotating with respect to the ICRS, so the rotation between them is the
identity. The two types are distinct because the *origin* differs and because
`GCRF` is where the Earth precession-nutation chain begins, not because the
axes are oriented differently.

See also: [`ICRF`](@ref).

# Example

```jldoctest
using AstroUniverse: earth
CoordinateSystem(earth, GCRF())

# output
CoordinateSystem:
  origin = Earth
  axes   = GCRF
```
"""
struct GCRF <: AbstractAxes end

"""
    CIRS()

Celestial Intermediate Reference System axes.

Earth-related IAU-2006 intermediate inertial frame using the Celestial
Intermediate Origin (CIO) convention. **Requires an Earth origin.**

# Example

```jldoctest
using AstroUniverse: earth
CoordinateSystem(earth, CIRS())

# output
CoordinateSystem:
  origin = Earth
  axes   = CIRS
```
"""
struct CIRS <: AbstractAxes end

"""
    TIRS()

Terrestrial Intermediate Reference System axes.

Earth-fixed IAU-2006 intermediate frame using the CIO convention. **Requires
an Earth origin.**

# Example

```jldoctest
using AstroUniverse: earth
CoordinateSystem(earth, TIRS())

# output
CoordinateSystem:
  origin = Earth
  axes   = TIRS
```
"""
struct TIRS <: AbstractAxes end

"""
    ITRF()

International Terrestrial Reference System axes.

Earth-fixed frame with polar-motion corrections from IERS Earth Orientation
Parameters. Production-grade Earth-fixed frame. **Requires an Earth origin.**

# Example

```jldoctest
using AstroUniverse: earth
CoordinateSystem(earth, ITRF())

# output
CoordinateSystem:
  origin = Earth
  axes   = ITRF
```
"""
struct ITRF <: AbstractAxes end

# --- Historical (FK5 / IAU-1976 family) -------------------------------------

"""
    MJ2000Eq()

Mean equator and equinox of MJ2000Eq (FK5) axes.

Origin-agnostic inertial. The axes describe a fixed spatial orientation
defined by Earth's mean equator and dynamical equinox at epoch MJ2000Eq.0;
they can legally be used with any origin.

# Example

```jldoctest
using AstroUniverse: earth
CoordinateSystem(earth, MJ2000Eq())

# output
CoordinateSystem:
  origin = Earth
  axes   = MJ2000Eq
```
"""
struct MJ2000Eq <: AbstractAxes end

"""
    MODEq()

Mean of Date, Equatorial (FK5) axes.

Earth-related inertial because the equator of date is Earth-specific. Requires
an Earth origin.

# Example

```jldoctest
using AstroUniverse: earth
CoordinateSystem(earth, MODEq())

# output
CoordinateSystem:
  origin = Earth
  axes   = MODEq
```
"""
struct MODEq <: AbstractAxes end

"""
    TODEq()

True of Date, Equatorial (FK5) axes.

Earth-related inertial because the equator of date is Earth-specific. Requires
an Earth origin.

# Example

```jldoctest
using AstroUniverse: earth
CoordinateSystem(earth, TODEq())

# output
CoordinateSystem:
  origin = Earth
  axes   = TODEq
```
"""
struct TODEq <: AbstractAxes end

"""
    MODEc()

Mean of Date, Ecliptic (FK5) axes.

Earth-orbit-related inertial because the ecliptic of date is Earth's orbital
plane. Requires an Earth origin.

# Example

```jldoctest
using AstroUniverse: earth
CoordinateSystem(earth, MODEc())

# output
CoordinateSystem:
  origin = Earth
  axes   = MODEc
```
"""
struct MODEc <: AbstractAxes end

"""
    TODEc()

True of Date, Ecliptic (FK5) axes.

Earth-orbit-related inertial because the ecliptic of date is Earth's orbital
plane. Requires an Earth origin.

# Example

```jldoctest
using AstroUniverse: earth
CoordinateSystem(earth, TODEc())

# output
CoordinateSystem:
  origin = Earth
  axes   = TODEc
```
"""
struct TODEc <: AbstractAxes end

"""
    PEF()

Pseudo Earth Fixed (FK5) axes.

Earth-fixed FK5 intermediate with Earth rotation applied through GAST and polar
motion omitted. `ITRF` includes the polar-motion correction. Requires an Earth
origin.

# Example

```jldoctest
using AstroUniverse: earth
CoordinateSystem(earth, PEF())

# output
CoordinateSystem:
  origin = Earth
  axes   = PEF
```
"""
struct PEF <: AbstractAxes end

"""
    TEME()

True Equator, Mean Equinox (FK5) axes.

The frame SGP4 propagates in, and therefore the frame a two-line element set
resolves to. It shares its equator with `TODEq` and differs from it by the
equation of the equinoxes: true equator with the mean equinox.

TEME supports TLE interoperability. Results should be converted to another
frame after propagation. Requires an Earth origin.

# Example

```jldoctest
using AstroUniverse: earth
CoordinateSystem(earth, TEME())

# output
CoordinateSystem:
  origin = Earth
  axes   = TEME
```
"""
struct TEME <: AbstractAxes end

# --- Static ecliptic ---------------------------------------------------------

"""
    MJ2000Ec()

Mean ecliptic and equinox of MJ2000Eq axes.

Origin-agnostic inertial. Static rotation from `MJ2000Eq` about the X-axis by
the mean obliquity of MJ2000Eq (ε₀ = 23.4392911° per IAU 1976/FK5).

# Example

```jldoctest
using AstroUniverse: earth
CoordinateSystem(earth, MJ2000Ec())

# output
CoordinateSystem:
  origin = Earth
  axes   = MJ2000Ec
```
"""
struct MJ2000Ec <: AbstractAxes end

# --- Body-fixed (Moon; SPICE-backed) ----------------------------------------

"""
    MoonPA()

The Moon principal-axes frame is used for dynamics.

The axes align with the Moon's principal axes of inertia. JPL lunar libration
integration produces this frame directly, and lunar gravity-field coefficients
are expressed in it.

Lunar maps, landing sites, surface features, and latitude-longitude coordinates
use [`MoonME`](@ref) instead. The frames differ by about 875 m at the surface.

Carries the full libration: physical, forced and free. Requires a Moon origin,
and requires `moon_pa_de440_200625.bpc` and `moon_de440_250416.tf` to be
loaded.

See also: [`MoonME`](@ref).

# Example

```jldoctest
using AstroUniverse: moon
CoordinateSystem(moon, MoonPA())

# output
CoordinateSystem:
  origin = Moon
  axes   = MoonPA
```
"""
struct MoonPA <: AbstractAxes end

"""
    MoonME()

The Moon mean Earth/mean rotation-axis frame is used for mapping.

The X axis points toward the mean direction of Earth, and the Z axis follows the
mean rotation axis. Published lunar coordinates, including landing sites,
craters, LOLA elevation grids, LRO imagery, and surface latitude and longitude,
use this frame.

For evaluating a lunar gravity field, use [`MoonPA`](@ref) instead. The two
differ by about 875 m on the surface, so the choice is not cosmetic.

The frame is separated from [`MoonPA`](@ref) by a fixed rotation of about 104
arcseconds and carries the same libration. It is not a smoothed version of the
principal axes; the two frames rotate together with different orientations.

Requires a Moon origin, and the same kernels as [`MoonPA`](@ref).

See also: [`MoonPA`](@ref).

# Example

```jldoctest
using AstroUniverse: moon
CoordinateSystem(moon, MoonME())

# output
CoordinateSystem:
  origin = Moon
  axes   = MoonME
```
"""
struct MoonME <: AbstractAxes end

# --- Body-fixed (Sun + planets; IAU 2015) -----------------------------------

"""
    CelestialBodyFixed{NAIFID} <: AbstractAxes

Body-fixed rotating axes for Sun, Mercury, Venus, Mars, Jupiter, Saturn,
Uranus, Neptune, or Pluto, using the IAU 2015 planet-rotation formulas from
AstroUniverse's `iau2015_orientation`. Type parameter `NAIFID::Int` records
the intended origin body; the origin-coupling rule requires the coordinate
system's origin to have a matching NAIF ID.

Not applicable to Earth (use IAU 2006 / FK5 Earth-frame edges) or Moon (use
`MoonPA` / `MoonME`).

# Constructors
- `CelestialBodyFixed(body)` — explicit form; encodes the body's NAIF ID in the type.
- `CelestialBodyFixed()` — unresolved sentinel; only valid when passed to
  `CoordinateSystem`, which fills in the NAIF ID from the origin.

# Examples
```julia
# Recommended: the origin supplies the body.
mars_fixed = CoordinateSystem(mars, CelestialBodyFixed())

# Equivalent, explicit:
mars_fixed = CoordinateSystem(mars, CelestialBodyFixed(mars))

# Direct axes-pair use (no CoordinateSystem) requires the explicit form:
M = axes_rotation(ICRF(), CelestialBodyFixed(mars), jd_tdb)
```
"""
struct CelestialBodyFixed{NAIFID} <: AbstractAxes end

CelestialBodyFixed(body::AbstractPoint) = CelestialBodyFixed{Int(body.naifid)}()
CelestialBodyFixed() = CelestialBodyFixed{0}()

"""
    naifid(::CelestialBodyFixed{N}) -> Int

Recover the NAIF ID stored in the axes type parameter.
"""
naifid(::CelestialBodyFixed{N}) where {N} = N

# NAIF ID → human-readable body name, for show + error messages.
const _BODY_NAME = Dict{Int,String}(
    10  => "Sun",     199 => "Mercury", 299 => "Venus",
    499 => "Mars",    599 => "Jupiter", 699 => "Saturn",
    799 => "Uranus",  899 => "Neptune", 999 => "Pluto",
)

_iau_body_name(n::Integer) = n == 0 ? "unresolved" : get(_BODY_NAME, Int(n), "NAIF-$(Int(n))")

function Base.show(io::IO, ::CelestialBodyFixed{N}) where {N}
    print(io, "CelestialBodyFixed(", _iau_body_name(N), ")")
end

Base.show(io::IO, ::MIME"text/plain", a::CelestialBodyFixed) = show(io, a)

# --- Orbit-relative ---------------------------------------------------------

# These three are one construction with different choices of primary direction,
# secondary direction, and labelling. They are kept as separate tags rather
# than one parametric tag because they are separate frames in every document a
# user will read, and because the sign conventions differ in ways a parameter
# would hide. `orbit_relative.jl` holds the shared construction.
#
# The set is deliberately closed at three. RSW and QSW are `RIC` under other
# names; TNW and NTW are `VNB` under other names. Providing aliases would mean
# four spellings of two frames.

"""
    RIC

Radial / in-track / cross-track axes of a reference orbit.

`R` follows the radius, `C` follows the orbit normal `r × v`, and `I` completes
the right-handed set. The `I` axis is exactly along-track only for a circular
orbit.

The standard frame for expressing relative position and orbit uncertainty.
Called RSW and QSW elsewhere; those are this frame, not different ones.

# Fields
None. The reference orbit passed to the transformation defines the frame.

# Notes
Needs a reference orbit, which AstroFrames cannot obtain on its own: pass it as
`reference_state`, a six-element inertial state in km and km/s. Acceleration is
optional here and, when given, includes the orbit plane turning under
out-of-plane force.

# Example
```julia
p = (; reference_state = [-4550.0, 2220.0, 4980.0, -3.10, -6.60, 0.12])

M = axes_rotation(ICRF(), RIC(), epoch, p)
separation_ric = M * separation_icrf      # relative state, radial/in-track/cross-track
```
"""
struct RIC <: AbstractAxes end

"""
    LVLH

Local-vertical / local-horizontal axes of a reference orbit.

AstroFrames uses the nadir-pointing convention: `z` is nadir (`−r̂`), `y` is the
negative orbit normal (`−ĥ`), and `x = y × z` completes the set and lies nearly
along-track. Other tools may use a different order or sign convention.

# Fields
None. See [`RIC`](@ref).

# Notes
Needs a reference orbit, passed as `reference_state`; acceleration is optional.

# Example
```julia
M = axes_rotation(ICRF(), LVLH(), epoch, (; reference_state = state))
M[3, :]    # the nadir direction in inertial axes
```
"""
struct LVLH <: AbstractAxes end

"""
    VNB

Velocity / normal / binormal axes of a reference orbit.

`V` along velocity, `N` along the orbit normal, `B = V × N`. The natural frame
for manoeuvres, whose components are usually stated along and across the
velocity.

# Fields
None. See [`RIC`](@ref).

# Notes
This frame turns for a different reason than [`RIC`](@ref) and [`LVLH`](@ref)
do. Those are built from the radius, whose direction swings around as the
spacecraft travels, so they turn with the orbit whatever forces are acting.
`V` follows the velocity, and its direction turns only under acceleration.

Acceleration is therefore a separate, optional input. Without
`reference_accel`, VNB describes an unaccelerated reference and has zero frame
rate. This form is appropriate for resolving a vector into instantaneous
along-track, normal, and binormal components. Supplying `reference_accel` in
km/s² includes the frame rate of an accelerating reference.

Called TNW and NTW elsewhere.

# Example
```julia
p = (; reference_state = state)

M = axes_rotation(ICRF(), VNB(), epoch, p)
Δv_inertial = M[1:3, 1:3]' * [0.1, 0.0, 0.0]    # 0.1 km/s along-track burn
```
"""
struct VNB <: AbstractAxes end

# --- Ambient-inertial marker -------------------------------------------------

"""
    Inertial

Marker for "the ambient inertial frame" — the inertial frame a state is
already expressed in, rather than any specific one.

Used where a quantity is defined relative to whatever inertial frame is in
play rather than to a named frame: an impulsive manoeuvre with
`axes=Inertial()` applies its ΔV components directly, with no rotation.
Contrast `VNB`, which rotates into an orbit-relative triad.

This is not a synonym for `ICRF` or `MJ2000Eq`. Naming a specific frame here
would be wrong whenever the state is expressed in a different inertial frame.

# Example

```jldoctest
needs_reference_orbit(Inertial())

# output
false
```
"""
struct Inertial <: AbstractAxes end


# =============================================================================
# Deprecated bindings
#
# The frames rename replaced these. `ICRFAxes` and `MJ2000Axes` were exported by
# the registered AstroFrames, so scripts written against them exist; they keep
# working for one release and warn on use rather than breaking on upgrade.
#
# Remove in the release after next.
# =============================================================================

Base.@deprecate_binding ICRFAxes   ICRF
Base.@deprecate_binding MJ2000Axes MJ2000Eq
