# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# State-transform dispatch API
#
# Public entry point:
#
#   axes_rotation(source_axes, target_axes, epoch[, params]) -> SMatrix{6,6,Float64,36}
#
# The 6×6 matrix `M` returned matches the SPICE `sxform` convention:
#   state_target = M * state_source            (state = [r; v])
#
# Each edge is a method on a pair of axes taking `EpochScales`. Routing between
# axes that share no edge composes edges along the route for the active frame
# theory. `epoch` arrives as a `Time` or a TDB Julian date and is converted to
# `EpochScales` once, at the public entry points at the end of this file.
# Static rotations ignore the epoch.
# =============================================================================

using SPICE: pxform as _spice_pxform, sxform as _spice_sxform, namfrm as _spice_namfrm
using SatelliteToolboxTransformations: DCM, r_mj2000_to_gcrf_iau2006, r_gcrf_to_mod_fk5, r_mod_to_tod_fk5, r_tod_to_pef_fk5, r_pef_to_itrf_fk5,
    r_gcrf_to_cirs_iau2006, r_tirs_to_itrf_iau2006,
    nutation_fk5, r_tod_to_teme,
    EARTH_ANGULAR_SPEED
using AstroEpochs: Time, TDB, JD
using AstroUniverse: AbstractFrameTheory, FK5, IAU2006, eop, frame_theory

const _J2000_TDB_JD = 2451545.0

_et_from_jd_tdb(jd_tdb::Real) = (jd_tdb - _J2000_TDB_JD) * 86400.0

# --- Epoch handling ---------------------------------------------------------
#
# `Time` is the public epoch type: it carries its own scale, so a caller
# cannot silently pass UTC where TDB was meant, and it keeps the two-part
# Julian date instead of crushing it into one Float64.
#
# Scales are derived **once**, at the boundary, and passed to the per-edge
# methods as plain numbers. That keeps `Time` — which is not `isbits` — out
# of the inner loop, and stops each edge from repeating the conversion.
#
# `ut1` is absent because no edge needs it yet; it arrives with the first
# Earth-spin edge, and requires EOP.

"""
    EpochScales

`EpochScales` holds an epoch in the time scales used by frame transformations.

Precession and nutation use terrestrial time, Earth rotation uses universal
time, and planetary orientation uses barycentric dynamical time. AstroFrames
derives the scales once for each transformation and passes them to each edge.

# Fields
Frame extensions read these values with [`epoch_tdb`](@ref), [`epoch_tt`](@ref),
and [`epoch_utc`](@ref) rather than accessing the fields directly.

- `tdb` — barycentric dynamical time, Julian date.
- `tt` — terrestrial time, Julian date.
- `utc` — coordinated universal time, Julian date.

# Notes
Built once per transform, at the boundary where a `Time` arrives, and passed
down as plain numbers. Deriving a scale costs more than most of the rotations
do, so an edge must not convert one itself.

# Example
```julia
function AstroFrames.axes_rotation(::ICRF, ::MyAxes, e::EpochScales)
    days = epoch_tdb(e) - 2451545.0
    ...
end
```
"""
struct EpochScales{T<:Real}
    tdb::T
    tt::T
    utc::T
end

@inline _scales(t::Time) = EpochScales(t.tdb.jd, t.tt.jd, t.utc.jd)

"""
    epoch_tdb(e::EpochScales)

Barycentric dynamical time, as a Julian date.

# Arguments
- `e::EpochScales` — the epoch bundle handed to a frame's `axes_rotation`.

# Returns
The TDB Julian date. This is the scale planetary orientation and ephemerides
are computed in, and the one a body-fixed frame wants.

# Notes
Frame models must use their specified time scale. TDB and TT differ by less
than 2 ms, while TT and UTC differ by more than a minute; substituting one for
another rotates a spinning frame by the scale error multiplied by its rate.

# Example
```julia
function AstroFrames.axes_rotation(::ICRF, ::MyAxes, e::EpochScales)
    days_since_j2000 = epoch_tdb(e) - 2451545.0
    ...
end
```
"""
@inline epoch_tdb(e::EpochScales) = e.tdb

"""
    epoch_tt(e::EpochScales)

Terrestrial time, as a Julian date.

# Arguments
- `e::EpochScales` — the epoch bundle handed to a frame's `axes_rotation`.

# Returns
The TT Julian date. This is the scale precession and nutation are computed in.

# Notes
TT runs 32.184 s ahead of atomic time and does not have leap seconds, which is
why the Earth-orientation series are written against it. See
[`epoch_tdb`](@ref) for the shape of a frame that uses one of these.

# Example

```julia
function AstroFrames.axes_rotation(::ICRF, ::MyAxes, e::EpochScales)
    jd = epoch_tt(e)
    ...
end
```
"""
@inline epoch_tt(e::EpochScales) = e.tt

"""
    epoch_utc(e::EpochScales)

Coordinated universal time, as a Julian date.

# Arguments
- `e::EpochScales` — the epoch bundle handed to a frame's `axes_rotation`.

# Returns
The UTC Julian date. This is the scale Earth orientation parameters are
tabulated against, so any frame reading polar motion or UT1−UTC wants it.

# Notes
UTC has leap seconds, so an interval computed by subtracting two of these is
not a physical duration. Use it to index the IERS tables, not to propagate.

See [`epoch_tdb`](@ref) for the shape of a frame that reads one of these.

The value is carried in `EpochScales` so each EOP-consuming edge uses the same
conversion.

# Example

```julia
function AstroFrames.axes_rotation(::ICRF, ::MyAxes, e::EpochScales)
    jd = epoch_utc(e)
    ...
end
```
"""
@inline epoch_utc(e::EpochScales) = e.utc

"""
    _ut1(e::EpochScales, jd_utc) -> Real

UT1 Julian date, from UTC plus the tabulated UT1−UTC.

UT1 is derived from the measured offset distributed by the IERS and stored in
the EOP tables. Earth rotates at 7.3e-5 rad/s, so one second of UT1 error is
7.3e-5 rad of longitude, or about 460 m at the equator.
"""
@inline _ut1(::EpochScales, jd_utc::Real) =
    jd_utc + eop(FK5()).Δut1_utc(jd_utc) / 86_400
# `float` rather than the value as given: an integer Julian date is a perfectly
# reasonable thing to write, and `Time` cannot represent one — it splits the
# date into two parts and the split is fractional. Without this,
# `axes_rotation(ICRF(), ITRF(), 2451545)` died with `InexactError: Int64(0.5)`
# several calls down, from a signature that says it takes any `Real`.
@inline _scales(jd_tdb::Real) = _scales(float(jd_tdb))
@inline _scales(jd_tdb::AbstractFloat) = _scales(Time(jd_tdb, zero(jd_tdb), :tdb, :jd))

# --- Identity ---------------------------------------------------------------

"""
    axes_rotation(source_axes, target_axes, epoch) -> SMatrix{6,6,Float64,36}
    axes_rotation(source_axes, target_axes, epoch, params) -> SMatrix{6,6,Float64,36}

The 6×6 transform that takes a Cartesian state expressed in `source_axes` to the
same state expressed in `target_axes`, as `state_target = M * state_source` with
the state `[r; v]` in km and km/s.

# Arguments
- `source_axes::AbstractAxes` — the axes the state is expressed in, such as `ICRF()`.
- `target_axes::AbstractAxes` — the axes to express it in, such as `ITRF()`.
- `epoch` — an `AstroEpochs.Time`, which carries its own time scale, or a Julian
  date in TDB.
- `params::NamedTuple` — evaluated numbers for axes defined by a reference orbit.
  [`RIC`](@ref), [`LVLH`](@ref) and [`VNB`](@ref) read `reference_state`, the
  reference orbit's position and velocity in ICRF, in km and km/s, and accept an
  optional `reference_accel` in km/s². Axes that need no parameters ignore them.

# Notes
The matrix has the block form of SPICE `sxform`: `M[1:3,1:3]` and `M[4:6,4:6]`
hold the rotation `R`, `M[4:6,1:3]` holds its time derivative `Ṙ`, and
`M[1:3,4:6]` is zero. For Earth frames `Ṙ` includes sidereal rotation and treats
precession, nutation, polar motion and obliquity of date as constant, as GMAT and
Vallado do; the omitted rates are about 0.8 mm/s at LEO and 4.6 mm/s at GEO. VNB
has zero rate when `reference_accel` is not given.

The route between two axes follows the active frame theory, IAU 2006 by default,
set with `set_frame_theory!`. A pair that can only be reached through the other
Earth theory is still transformed, and warns once per pair.

Axes defined by a reference orbit raise an `ArgumentError` naming the missing
field when `params` does not carry `reference_state`. A pair with no route between
them raises an `ArgumentError` listing the supported edges.

# Returns
The 6×6 transform as an `SMatrix{6,6,Float64,36}`.

# Example
```jldoctest
using AstroEpochs
epoch = Time("2020-01-01T00:00:00.000", UTC(), ISOT())
M = axes_rotation(ICRF(), ITRF(), epoch)

# The same state in Earth-fixed axes, and back again
x_icrf = [7000.0, 0.0, 0.0, 0.0, 7.546, 0.0]
x_itrf = M * x_icrf
round.(axes_rotation(ITRF(), ICRF(), epoch) * x_itrf; digits = 9) == x_icrf

# output
true
```
"""
function axes_rotation end

@inline function axes_rotation(::A, ::A, ::EpochScales) where {A<:AbstractAxes}
    return SMatrix{6,6,Float64,36}(I)
end

# --- ICRF ↔ CelestialBodyFixed{NAIFID} --------------------------------------
#
# How a body is oriented is a property of the body, held in `AstroUniverse`,
# so this edge asks rather than deciding. That is what makes the frame usable
# for a body this package has never heard of: register a model, and the edge
# below works unchanged.
#
# `NAIFID` is baked into the axes type, so the lookup is a type-level constant
# and the model resolves through one function barrier.
#
# For the reverse: the inverse of a state transform of block form `[R 0; Ṙ R]`
# is `[Rᵀ 0; Ṙᵀ Rᵀ]` (verifiable by direct block multiplication — the
# off-diagonal cross-terms cancel because R is orthogonal and ṘᵀR + RᵀṘ =
# d/dt(RᵀR) = d/dt(I) = 0).

function axes_rotation(::ICRF, ::CelestialBodyFixed{N}, e::EpochScales) where {N}
    return _body_fixed_rotation(orientation_model(N), N, e)
end

# `orientation_model` returns an `AbstractOrientationModel`, because which model
# a body uses is a runtime setting. The barrier confines that to one call: `M`
# is concrete inside, so the model's own method resolves statically and its
# arithmetic inlines. Same shape as `_routed_rotation` does for the frame
# theory, and for the same reason.
#
# The dispatch itself is not free — see the note on the method above.
@inline _body_fixed_rotation(model::M, naifid, e::EpochScales) where {M<:AbstractOrientationModel} =
    body_axes_rotation(model, naifid, epoch_tdb(e))

function axes_rotation(::CelestialBodyFixed{N}, ::ICRF, e::EpochScales) where {N}
    return _invert_rotation(axes_rotation(ICRF(), CelestialBodyFixed{N}(), e))
end

# --- ICRF ↔ MoonPA / MoonME (SPICE-backed) ----------------------------------
#
# NOTE: the SPICE frame named "J2000" is ICRF-aligned — it is NOT our
# `MJ2000Eq` (FK5 mean equator and equinox), which differs from it by the
# ~23 mas frame bias. Mapping SPICE "J2000" to `ICRF` here is deliberate and
# correct. Do not "fix" these string literals to match our tag names.
#
# These need lunar kernels that are not among the defaults, so the common way
# to meet them is without the kernels loaded. Left alone, SPICE reports "The
# frame MOON_ME was not recognized as a known reference frame" — naming a
# string the user never typed, from a package they did not call, with nothing
# about what to load. `_lunar_sxform` turns that into an error that says.

"""
    _lunar_frame_error(frame, et, cause) -> ArgumentError

Why a lunar frame lookup failed, in terms of what to do about it.

`namfrm` returns 0 for a frame no loaded kernel defines, which separates the
two causes: no frame kernel at all, or a frame kernel with no orientation data
covering this epoch. They need different fixes, so they get different messages.
"""
function _lunar_frame_error(frame::AbstractString, et::Real, cause)
    axes = frame == "MOON_PA" ? "MoonPA" : "MoonME"
    load = """
               download_spice_kernel("moon_pa_de440_200625.bpc")
               download_spice_kernel("moon_de440_250416.tf")
               load_spice_kernel("moon_pa_de440_200625.bpc")
               load_spice_kernel("moon_de440_250416.tf")"""

    if _spice_namfrm(String(frame)) == 0
        return ArgumentError(
            "$(axes) axes need the lunar frame kernels, and none is loaded.
" *
            "
These are loaded for you at startup, so reaching this usually " *
            "means they were unloaded — `unload_all_spice_kernels()` clears " *
            "them. Reload with:

" * load *
            "

The frame kernel must be moon_de440_250416.tf. The older " *
            "moon_080317.tf names the DE421 frame and will not work against " *
            "DE440 orientation data.
" *
            "
SPICE reported: $(sprint(showerror, cause))")
    end

    return ArgumentError(
        "$(axes) axes are defined, but no lunar orientation data covers " *
        "$(Time(et / 86_400 + 2451545.0, 0.0, TDB(), JD()).isot) TDB.
" *
        "
moon_pa_de440_200625.bpc covers 1550 through 2650. If your epoch " *
        "is inside that, the kernel is probably not loaded:

" * load *
        "

SPICE reported: $(sprint(showerror, cause))")
end

"""`sxform` for the lunar frames, reporting a missing kernel as such."""
function _lunar_sxform(from::AbstractString, to::AbstractString, et::Real)
    try
        return _spice_sxform(String(from), String(to), Float64(et))
    catch cause
        # Which of the two is the lunar frame — the other end is always J2000.
        lunar = from == "J2000" ? to : from
        throw(_lunar_frame_error(lunar, et, cause))
    end
end

function axes_rotation(::ICRF, ::MoonPA, e::EpochScales)
    et = _et_from_jd_tdb(e.tdb)
    M = _lunar_sxform("J2000", "MOON_PA", et)
    return SMatrix{6,6,Float64,36}(M)
end

function axes_rotation(::MoonPA, ::ICRF, e::EpochScales)
    et = _et_from_jd_tdb(e.tdb)
    M = _lunar_sxform("MOON_PA", "J2000", et)
    return SMatrix{6,6,Float64,36}(M)
end

function axes_rotation(::ICRF, ::MoonME, e::EpochScales)
    et = _et_from_jd_tdb(e.tdb)
    M = _lunar_sxform("J2000", "MOON_ME", et)
    return SMatrix{6,6,Float64,36}(M)
end

function axes_rotation(::MoonME, ::ICRF, e::EpochScales)
    et = _et_from_jd_tdb(e.tdb)
    M = _lunar_sxform("MOON_ME", "J2000", et)
    return SMatrix{6,6,Float64,36}(M)
end

# --- MJ2000Eq ↔ MODEq (IAU-1976 precession) ----------------------------------
#
# STB names the inertial end of its FK5 precession "GCRF", but its own
# docstring is explicit that without EOP corrections this is "what is usually
# called the J2000 reference frame" — our `MJ2000Eq`, not our `GCRF`.
#
# Needs TT only: no EOP, no UT1. `Ṙ` is neglected (see `_rotation_no_rate`).

"""
    axes_rotation(::MJ2000Eq, ::MODEq, epoch) -> SMatrix{6,6,Float64,36}

IAU-1976 precession from the FK5 mean equator and equinox of J2000 to the
mean equator and equinox of date.

Requires no Earth orientation data. `epoch` may be a `Time` or a TDB Julian
date; the TT needed by the model is derived internally.
"""
axes_rotation(::MJ2000Eq, ::MODEq, e::EpochScales) =
    _rotation_no_rate(r_gcrf_to_mod_fk5(DCM, e.tt))

axes_rotation(::MODEq, ::MJ2000Eq, e::EpochScales) =
    _invert_rotation(axes_rotation(MJ2000Eq(), MODEq(), e))

# --- Which theory an edge belongs to ----------------------------------------
#
# Each edge declares its own theory, and the theory selects the EOP table
# (`AstroUniverse.eop`). This is what makes passing the wrong table
# unreachable rather than merely discouraged: an edge cannot ask for a table
# by name, only for its own theory's.
#
# `nothing` means the edge needs no Earth orientation data at all.
#
# A route whose edges span both concrete theories is a modelling choice, not
# an error — it warns once (`FR-FRAME-4`) rather than failing.

"""
    edge_theory(source_axes, target_axes) -> AbstractFrameTheory | Nothing

Earth precession-nutation theory an edge belongs to, or `nothing` if it needs
no Earth orientation data.

# Arguments
- `source_axes::AbstractAxes` — one end of the edge.
- `target_axes::AbstractAxes` — the other end. The theory is the same in both
  directions.

# Returns
The theory this edge is defined under — `FK5()` or `IAU2006()` — or `nothing`
where the edge reads no Earth orientation data at all. The theory selects the
EOP table, which is what makes passing the wrong one unreachable rather than
merely discouraged.

# Example

```jldoctest
using AstroUniverse: FK5, IAU2006
edge_theory(MJ2000Eq(), MODEq()), edge_theory(GCRF(), CIRS()), edge_theory(ICRF(), MJ2000Eq())

# output
(FK5(), IAU2006(), nothing)
```
"""
edge_theory(::AbstractAxes, ::AbstractAxes) = nothing

edge_theory(::MJ2000Eq, ::MODEq) = FK5()
edge_theory(::MODEq, ::MJ2000Eq) = FK5()
edge_theory(::MODEq, ::TODEq)    = FK5()
edge_theory(::TODEq, ::MODEq)    = FK5()
edge_theory(::TODEq, ::PEF)      = FK5()
edge_theory(::PEF, ::TODEq)      = FK5()
edge_theory(::PEF, ::ITRF)       = FK5()
edge_theory(::ITRF, ::PEF)       = FK5()
edge_theory(::MODEq, ::MODEc)    = FK5()
edge_theory(::MODEc, ::MODEq)    = FK5()
edge_theory(::TODEq, ::TODEc)    = FK5()
edge_theory(::TODEc, ::TODEq)    = FK5()
edge_theory(::TODEq, ::TEME)     = FK5()
edge_theory(::TEME, ::TODEq)     = FK5()

edge_theory(::GCRF, ::CIRS)      = IAU2006()
edge_theory(::CIRS, ::GCRF)      = IAU2006()
edge_theory(::CIRS, ::TIRS)      = IAU2006()
edge_theory(::TIRS, ::CIRS)      = IAU2006()
edge_theory(::TIRS, ::ITRF)      = IAU2006()
edge_theory(::ITRF, ::TIRS)      = IAU2006()

# --- EOP unit conversion ----------------------------------------------------
#
# IERS distributes polar motion in arcseconds and the celestial-pole
# corrections in milliarcseconds; the rotation models want radians. STB's
# high-level wrappers apply these internally, but we call the per-edge `r_*`
# functions directly, so the conversion is ours to do.

const _ARCSEC_TO_RAD      = deg2rad(1 / 3600)
const _MILLIARCSEC_TO_RAD = deg2rad(1 / 3_600_000)

# --- PEF ↔ ITRF (polar motion) -----------------------------------------------
#
# Takes no epoch of its own: `r_pef_to_itrf_fk5` consumes the pole coordinates
# directly, and all the time dependence sits in the EOP lookup.
#
# `Ṙ` is neglected. Polar motion wanders by a few tenths of an arcsecond over
# ~430 days (the Chandler wobble plus an annual term), so its rate is ~1e-11
# rad/s — the same order as precession, and seven orders below the Earth spin
# that dominates the neighbouring edge.

"""
    axes_rotation(::PEF, ::ITRF, epoch) -> SMatrix{6,6,Float64,36}

Polar motion from the pseudo Earth-fixed frame to the terrestrial reference
frame, using the IERS pole coordinates.

This is the edge that distinguishes `PEF` from `ITRF`: `PEF` has Earth
rotation applied and polar motion not.

`epoch` may be a `Time` or a TDB Julian date.
"""
function axes_rotation(::PEF, ::ITRF, e::EpochScales)
    table  = eop(FK5())
    jd_utc = epoch_utc(e)
    # IERS distributes pole coordinates in arcseconds.
    x_p = table.x(jd_utc) * _ARCSEC_TO_RAD
    y_p = table.y(jd_utc) * _ARCSEC_TO_RAD
    return _rotation_no_rate(r_pef_to_itrf_fk5(DCM, x_p, y_p))
end

function axes_rotation(::ITRF, ::PEF, e::EpochScales)
    return _invert_rotation(axes_rotation(PEF(), ITRF(), e))
end

# --- TODEq ↔ TEME (equation of the equinoxes) --------------------------------
#
# TEME shares its equator with TODEq and differs only in the equinox: true
# equator, *mean* equinox. The rotation between them is the equation of the
# equinoxes, so there is no spin and `Ṙ` is zero.
#
# Attaching TEME here rather than to PEF is deliberate. STB also offers a
# direct TEME ↔ PEF rotation built on GMST, and since GAST = GMST + the
# equation of the equinoxes, routing TEME → TODEq → PEF must agree with it.
# That equivalence is a free cross-check, and it is tested.
#
# TEME exists for TLE interoperability. It is where SGP4 leaves a state, not a
# frame to work in.

"""
    axes_rotation(::TODEq, ::TEME, epoch) -> SMatrix{6,6,Float64,36}

True equator and equinox of date to the true-equator, mean-equinox frame, via
the equation of the equinoxes.

Applies the same IERS celestial-pole corrections as the nutation edge, so the
`TODEq` end matches the one that edge produces.
"""
function axes_rotation(::TODEq, ::TEME, e::EpochScales)
    table  = eop(FK5())
    jd_utc = epoch_utc(e)
    δΔϵ = table.δΔϵ(jd_utc) * _MILLIARCSEC_TO_RAD
    δΔψ = table.δΔψ(jd_utc) * _MILLIARCSEC_TO_RAD
    return _rotation_no_rate(r_tod_to_teme(DCM, e.tt, δΔϵ, δΔψ))
end

axes_rotation(::TEME, ::TODEq, e::EpochScales) =
    _invert_rotation(axes_rotation(TODEq(), TEME(), e))

# --- MODEq ↔ MODEc and TODEq ↔ TODEc (obliquity of date) ---------------------
#
# Equatorial to ecliptic is a single rotation about X by the obliquity. Which
# obliquity is the whole content of these two edges:
#
#   MODEc uses the MEAN obliquity ε_A — mean equator, mean ecliptic.
#   TODEc uses the TRUE obliquity ε_A + Δε — true equator, true ecliptic.
#
# The `δΔϵ` term matters for consistency, not precision. Our `MODEq ↔ TODEq`
# edge applies the IERS celestial-pole correction, so the TODEq this rotates
# out of is the *corrected* true-of-date frame. Using ε_A + Δε alone would
# rotate to an ecliptic belonging to a slightly different equator — an
# inconsistency of a few tens of milliarcseconds between two edges that are
# supposed to share a frame.
#
# `Ṙ` is neglected: the obliquity drifts ~47″/century, about 7e-11 rad/s.

# The mean obliquity is written out here rather than taken from
# `nutation_fk5`, which returns it from Vallado's decimal-degree form of the
# IAU-1980 expression:
#
#     ε = 23.439291° − 0.0130042 T − 1.64e−7 T² + 5.04e−7 T³
#
# The IERS conventions publish the same expression in arcseconds, and the two
# roundings are not the same number: 23.439291° is 84381.4476″, which is 0.4 mas
# below 84381.448″, and the linear terms differ by another 0.12 mas/century.
# Measured against ERFA the offset was 0.400 mas at J2000 drifting to 0.430 mas
# by 2025 — three orders of magnitude above every other edge.
#
# It shows up only here. In the nutation matrix an error in ε enters twice with
# opposite signs and cancels to first order, which is why `MODEq ↔ TODEq` agreed
# to 0.03 µas while the ecliptic edges, a single rotation by ε, took the whole
# offset.
#
# Δε still comes from `nutation_fk5`: the nutation in obliquity is a separate
# quantity from the mean obliquity, and that one agrees.

@inline function _mean_obliquity(e::EpochScales)
    T = (epoch_tt(e) - _J2000_TDB_JD) / 36_525
    return _ARCSEC_TO_RAD * @evalpoly(T, 84_381.448, -46.8150, -0.00059, 0.001813)
end

@inline function _true_obliquity(e::EpochScales)
    _, Δε, _ = nutation_fk5(epoch_tt(e))
    δΔϵ = eop(FK5()).δΔϵ(epoch_utc(e)) * _MILLIARCSEC_TO_RAD
    return _mean_obliquity(e) + Δε + δΔϵ
end

"""
    _earth_rotation_angle(jd_ut1::Real) -> Real

Earth Rotation Angle at `jd_ut1`, in radians, reduced to `[0, 2π)`.

# Notes
Written out here rather than taken from `SatelliteToolboxTransformations`,
because of how the reduction has to be arranged. ERA is

    θ = 2π (0.7790572732640 + 1.00273781191135448 d)

with `d` the days elapsed since J2000. Formed directly, the bracket reaches
about 9200 by 2025, where a `Float64` spacing is already 2.4 µas of angle. The
subsequent reduction then discards the whole turns after they have reduced the
precision. This arrangement lost 0.7 µas against ERFA by 2025, with the error
growing with epoch.

Splitting the date keeps the turns out of the multiply. Only the *excess* over
one turn per day, 0.00273781191135448, multiplies the elapsed days; the turns
themselves arrive through the fractional part of the Julian date, already in
`[0, 1)`. The argument stays O(20) rather than O(9200) and the result matches
ERFA to 0.001 µas across the IERS record.

This is IAU SOFA's own arrangement, in `eraEra00`.
"""
@inline function _earth_rotation_angle(jd_ut1::Real)
    # Split at the half day, so the fractional part is exact for the `.5`
    # Julian dates that dominate in practice.
    jd_hi = floor(jd_ut1 - 0.5) + 0.5
    jd_lo = jd_ut1 - jd_hi

    d = jd_ut1 - _J2000_TDB_JD             # elapsed days; only ever times 0.0027…
    f = mod(jd_hi, one(jd_hi)) + jd_lo     # fractional day, in [0, 1)

    return mod2pi(2 * pi * (f + 0.7790572732640 + 0.00273781191135448 * d))
end

"""
    axes_rotation(::MODEq, ::MODEc, epoch) -> SMatrix{6,6,Float64,36}

Transforms mean equator-of-date axes to mean ecliptic-of-date axes by rotating
about X through the IAU 1980 mean obliquity ε_A.
"""
axes_rotation(::MODEq, ::MODEc, e::EpochScales) =
    _rotation_no_rate(_Rx(_mean_obliquity(e)))

axes_rotation(::MODEc, ::MODEq, e::EpochScales) =
    _invert_rotation(axes_rotation(MODEq(), MODEc(), e))

"""
    axes_rotation(::TODEq, ::TODEc, epoch) -> SMatrix{6,6,Float64,36}

Transforms true equator-of-date axes to true ecliptic-of-date axes by rotating
about X through the true obliquity, `ε_A + Δε`. The rotation includes the IERS
`δΔϵ` correction used by the `TODEq` nutation edge.
"""
axes_rotation(::TODEq, ::TODEc, e::EpochScales) =
    _rotation_no_rate(_Rx(_true_obliquity(e)))

axes_rotation(::TODEc, ::TODEq, e::EpochScales) =
    _invert_rotation(axes_rotation(TODEq(), TODEc(), e))

# --- TODEq ↔ PEF (Earth rotation via GAST) -----------------------------------
#
# The first edge whose rate is *not* neglected. Earth spin dominates every
# other rate in the chain by seven orders of magnitude — 7.3e-5 rad/s against
# ~1e-11 for precession and nutation — so `Ṙ` here is the whole story for
# velocity, and `_rotation_with_spin` is used rather than `_rotation_no_rate`.
#
# Needs UT1 (Earth rotation angle), TT (the nutation the GAST expression
# depends on), and the FK5 EOP series for both Δut1_utc and the length-of-day
# correction to the spin rate.
#
# Polar motion is NOT applied: that is the PEF ↔ ITRF edge. PEF is precisely
# the frame with Earth rotation applied and polar motion not.

"""
    axes_rotation(::TODEq, ::PEF, epoch) -> SMatrix{6,6,Float64,36}

Earth rotation from the true equator and equinox of date to the pseudo
Earth-fixed frame, via Greenwich apparent sidereal time.

Unlike the precession and nutation edges, this one carries a real rotation
rate: `Ṙ = −[ω]ₓR` with `ω` the LOD-corrected Earth spin rate. Velocity
transforms accordingly.

This transformation does not apply polar motion; the `PEF ↔ ITRF` edge does.

`epoch` may be a `Time` or a TDB Julian date.
"""
function axes_rotation(::TODEq, ::PEF, e::EpochScales)
    table  = eop(FK5())
    jd_utc = epoch_utc(e)
    jd_ut1 = _ut1(e, jd_utc)
    δΔψ    = table.δΔψ(jd_utc) * _MILLIARCSEC_TO_RAD

    # LOD is distributed in milliseconds; a longer day is a slower Earth.
    ω = EARTH_ANGULAR_SPEED * (1 - table.lod(jd_utc) / 86_400_000)

    return _rotation_with_spin(r_tod_to_pef_fk5(DCM, jd_ut1, e.tt, δΔψ), ω)
end

function axes_rotation(::PEF, ::TODEq, e::EpochScales)
    return _invert_rotation(axes_rotation(TODEq(), PEF(), e))
end

# --- ICRF ↔ GCRF (identity) ---------------------------------------------------
#
# Not an approximation. The GCRS is *defined* as kinematically non-rotating
# with respect to the ICRS — obtained by a Lorentz transformation with no
# spatial rotation — so their axes are aligned and there is nothing to model.
# What separates them is the origin and the relativistic treatment of
# observations, neither of which is an axes rotation.
#
# Geodesic precession does not belong here. It concerns the kinematically
# non-rotating GCRS versus a *dynamically* non-rotating geocentric frame,
# which is a different comparison.

"""
    axes_rotation(::ICRF, ::GCRF, epoch) -> SMatrix{6,6,Float64,36}

Returns the identity transformation because `GCRF` shares its orientation with
`ICRF` by definition. The frames differ in origin, not in axes.
"""
axes_rotation(::ICRF, ::GCRF, ::EpochScales) = SMatrix{6,6,Float64,36}(I)
axes_rotation(::GCRF, ::ICRF, ::EpochScales) = SMatrix{6,6,Float64,36}(I)

# --- GCRF ↔ CIRS (IAU-2006 precession-nutation, CIO-based) -------------------
#
# Needs TT and the CIP offsets δx, δy from the IAU-2006 EOP series. `Ṙ` is
# neglected, as for the FK5 precession and nutation edges.

"""
    axes_rotation(::GCRF, ::CIRS, epoch) -> SMatrix{6,6,Float64,36}

IAU-2006/2010 precession-nutation from the GCRF to the Celestial Intermediate
Reference System, using the CIO-based formulation.

Applies the IERS CIP offsets from the IAU-2006 EOP series.
"""
function axes_rotation(::GCRF, ::CIRS, e::EpochScales)
    table  = eop(IAU2006())
    jd_utc = epoch_utc(e)
    δx = table.δx(jd_utc) * _MILLIARCSEC_TO_RAD
    δy = table.δy(jd_utc) * _MILLIARCSEC_TO_RAD
    return _rotation_no_rate(r_gcrf_to_cirs_iau2006(DCM, e.tt, δx, δy))
end

axes_rotation(::CIRS, ::GCRF, e::EpochScales) =
    _invert_rotation(axes_rotation(GCRF(), CIRS(), e))

# --- CIRS ↔ TIRS (Earth rotation via ERA) ------------------------------------
#
# The IAU-2006 counterpart of the FK5 GAST edge, and like it a genuine spin
# edge: `Ṙ` is the dominant velocity term, so `_rotation_with_spin` is used.
#
# Needs UT1 only — the Earth Rotation Angle is a linear function of UT1, which
# is what makes the CIO formulation simpler than GAST.

"""
    axes_rotation(::CIRS, ::TIRS, epoch) -> SMatrix{6,6,Float64,36}

Earth rotation from the Celestial to the Terrestrial Intermediate Reference
System, via the Earth Rotation Angle.

Carries the Earth rotation rate, like its FK5 counterpart `TODEq ↔ PEF`. This
transformation does not apply polar motion; the `TIRS ↔ ITRF` edge does.
"""
function axes_rotation(::CIRS, ::TIRS, e::EpochScales)
    table  = eop(IAU2006())
    jd_utc = epoch_utc(e)
    jd_ut1 = jd_utc + table.Δut1_utc(jd_utc) / 86_400
    ω = EARTH_ANGULAR_SPEED * (1 - table.lod(jd_utc) / 86_400_000)
    return _rotation_with_spin(_Rz(_earth_rotation_angle(jd_ut1)), ω)
end

axes_rotation(::TIRS, ::CIRS, e::EpochScales) =
    _invert_rotation(axes_rotation(CIRS(), TIRS(), e))

# --- TIRS ↔ ITRF (polar motion) ----------------------------------------------
#
# The IAU-2006 counterpart of PEF ↔ ITRF. Unlike the FK5 version this one does
# take an epoch: the IAU-2006 polar-motion matrix includes the TIO locator s′,
# which is a function of time. `Ṙ` is neglected.

"""
    axes_rotation(::TIRS, ::ITRF, epoch) -> SMatrix{6,6,Float64,36}

Polar motion from the Terrestrial Intermediate Reference System to the
terrestrial reference frame, IAU-2006 formulation.

This is the edge that distinguishes `TIRS` from `ITRF`: `TIRS` has Earth
rotation applied and polar motion not.
"""
function axes_rotation(::TIRS, ::ITRF, e::EpochScales)
    table  = eop(IAU2006())
    jd_utc = epoch_utc(e)
    x_p = table.x(jd_utc) * _ARCSEC_TO_RAD
    y_p = table.y(jd_utc) * _ARCSEC_TO_RAD
    return _rotation_no_rate(r_tirs_to_itrf_iau2006(DCM, e.tt, x_p, y_p))
end

axes_rotation(::ITRF, ::TIRS, e::EpochScales) =
    _invert_rotation(axes_rotation(TIRS(), ITRF(), e))

# --- MODEq ↔ TODEq (IAU-1980 nutation) ---------------------------------------
#
# Needs TT, plus the IERS celestial-pole corrections (δΔψ, δΔϵ) from the FK5
# EOP series. Those corrections are not cosmetic: they correct known errors in
# the 1980 nutation model, and they are larger than "small" suggests — δΔψ was
# about −108 mas at 2020, shifting the total nutation rotation by ~0.11″.
#
# STB's argument order is (jd_tt, δΔϵ, δΔψ) — obliquity before longitude,
# which is the reverse of how they are usually written. Swapping them is a
# silent error of the right magnitude, so they are passed by name below.
#
# `Ṙ` is neglected (see `_rotation_no_rate`).

"""
    axes_rotation(::MODEq, ::TODEq, epoch) -> SMatrix{6,6,Float64,36}

IAU-1980 nutation from the mean equator and equinox of date to the true
equator and equinox of date.

Applies the IERS celestial-pole corrections from the FK5 EOP series, so the
result is the corrected IAU-80 nutation rather than the original model.

`epoch` may be a `Time` or a TDB Julian date.
"""
function axes_rotation(::MODEq, ::TODEq, e::EpochScales)
    table = eop(FK5())
    jd_utc = epoch_utc(e)
    δΔϵ = table.δΔϵ(jd_utc) * _MILLIARCSEC_TO_RAD
    δΔψ = table.δΔψ(jd_utc) * _MILLIARCSEC_TO_RAD
    return _rotation_no_rate(r_mod_to_tod_fk5(DCM, e.tt, δΔϵ, δΔψ))
end

axes_rotation(::TODEq, ::MODEq, e::EpochScales) =
    _invert_rotation(axes_rotation(MODEq(), TODEq(), e))

# --- ICRF ↔ MJ2000Eq (frame bias, static) ------------------------------------
#
# The rotation between the ICRF and the FK5 mean equator and equinox of J2000.
# It is a *bias* — epoch-independent — of about 23 mas, which is small but not
# negligible: applied over an interplanetary baseline it is ~16 km at 1 AU.
#
# Taken from SatelliteToolboxTransformations rather than hand-derived. The
# bias is three small angles whose sign conventions differ between references,
# and a sign error here is invisible in magnitude while being wrong in
# direction. STB's implementation follows Vallado (2013).
#
# NOTE on names: STB calls the inertial end "GCRF" and the other end "MJ2000".
# STB's GCRF here is the ICRF-aligned top-level inertial frame — our `ICRF` —
# not our `GCRF`, which is Earth-restricted and carries geodesic precession.

const _M_ICRF_TO_MJ2000EQ = _rotation_no_rate(
    inv(r_mj2000_to_gcrf_iau2006(DCM, 0.0))
)

"""
    axes_rotation(::ICRF, ::MJ2000Eq, epoch) -> SMatrix{6,6,Float64,36}

Frame bias between the ICRF and the FK5 mean equator and equinox of J2000.

Epoch-independent: the argument is accepted for interface uniformity and
ignored. `Ṙ` is exactly zero.
"""
axes_rotation(::ICRF, ::MJ2000Eq, ::EpochScales) = _M_ICRF_TO_MJ2000EQ
axes_rotation(::MJ2000Eq, ::ICRF, ::EpochScales) = _invert_rotation(_M_ICRF_TO_MJ2000EQ)

# --- MJ2000Eq ↔ MJ2000Ec (static rotation about X by ε₀) ----------------------

# Mean obliquity of J2000 (the epoch), IAU 1976/FK5. 23°26'21.448″ = 23.4392911111°.
const _OBLIQUITY_J2000_RAD = deg2rad(23.4392911111)

const _M_J2000_TO_MJ2000EC = _rotation_no_rate(_Rx(_OBLIQUITY_J2000_RAD))

axes_rotation(::MJ2000Eq, ::MJ2000Ec, ::EpochScales) = _M_J2000_TO_MJ2000EC
axes_rotation(::MJ2000Ec, ::MJ2000Eq, ::EpochScales) = _invert_rotation(_M_J2000_TO_MJ2000EC)

# --- Unresolved CelestialBodyFixed sentinel --------------------------------------

_unresolved_body_msg() =
    "Unresolved `CelestialBodyFixed()` cannot be used directly in `axes_rotation`. " *
    "Pass an explicit body (e.g. `CelestialBodyFixed(mars)`), or construct via " *
    "`CoordinateSystem(mars, CelestialBodyFixed())` which fills in the body from the origin."

axes_rotation(::CelestialBodyFixed{0}, ::AbstractAxes, ::EpochScales) = throw(ArgumentError(_unresolved_body_msg()))
axes_rotation(::AbstractAxes, ::CelestialBodyFixed{0}, ::EpochScales) = throw(ArgumentError(_unresolved_body_msg()))
# Disambiguations: the ICRF↔CelestialBodyFixed{N} methods above match N=0 too.
axes_rotation(::ICRF, ::CelestialBodyFixed{0}, ::EpochScales)          = throw(ArgumentError(_unresolved_body_msg()))
axes_rotation(::CelestialBodyFixed{0}, ::ICRF, ::EpochScales)          = throw(ArgumentError(_unresolved_body_msg()))
axes_rotation(::CelestialBodyFixed{0}, ::CelestialBodyFixed{0}, ::EpochScales) = throw(ArgumentError(_unresolved_body_msg()))

# --- Fallback ---------------------------------------------------------------

# --- Routing -----------------------------------------------------------------
#
# Most pairs have no direct edge; they are reached by composing the chain. The
# intermediates are *declared*, not searched for at run time, so composition
# resolves at compile time and a routed pair costs what a hand-written one
# costs.
#
# The route depends on the frame theory for pairs whose endpoints are
# theory-neutral: ICRF and ITRF both exist in either branch, so the theory
# chooses which chain connects them. That is the theory setting's real job.

"""
    _route(source_axes, target_axes, theory) -> Tuple | Nothing

Intermediate axes connecting two frames with no direct edge, in order, or
`nothing` if no route is declared.
"""
_route(::AbstractAxes, ::AbstractAxes, ::AbstractFrameTheory) = nothing

"""
    _route_any(source_axes, target_axes) -> Tuple | Nothing

Returns a route without restricting the search to the active theory. AstroFrames
uses this fallback only when the active theory cannot reach the target and emits
the `FR-FRAME-4` warning.
"""
_route_any(::AbstractAxes, ::AbstractAxes) = nothing

# Warn-once bookkeeping. A transform inside a propagation loop must not warn on
# every step, so the key is the frame pair and the theory in force.
const _OUT_OF_THEORY_WARNED = Set{Tuple{DataType,DataType,DataType}}()
const _OUT_OF_THEORY_LOCK   = ReentrantLock()

function _warn_out_of_theory(source, target, theory)
    key = (typeof(source), typeof(target), typeof(theory))
    fresh = lock(_OUT_OF_THEORY_LOCK) do
        key in _OUT_OF_THEORY_WARNED ? false : (push!(_OUT_OF_THEORY_WARNED, key); true)
    end
    fresh || return nothing
    @warn string(
        nameof(typeof(source)), " → ", nameof(typeof(target)),
        " cannot be routed within the active frame theory ", nameof(typeof(theory)),
        ", so it is routed through edges of the other Earth-rotation theory. ",
        "The result composes links from more than one model. ",
        "Set the theory to match the frames you name (`set_frame_theory!`), or name frames ",
        "from the active theory's chain. Warned once per frame pair.")
    return nothing
end

# The graph is declared as **edges**, and the routes are searched for once, at
# load time. Declaring routes instead would be O(n²) statements that can drift
# out of step with the edges that actually exist — and would silently miss any
# pair the author did not think of, which is exactly what happened when this
# was first written as declared chains.
#
# The search runs here, not on the call path. That was the whole of the
# objection to the prototype's design: a BFS with a `Dict` lookup per edge at
# transform time is a dynamic dispatch on a hot path. Running the same search
# once and emitting fixed tuples keeps the generality and loses the cost.
#
# Each theory sees its own subgraph: theory-neutral edges are always available,
# and a theory's own edges are added to them. A pair reachable only by crossing
# theories therefore gets no route under either — which is correct, and is the
# `FR-FRAME-4` case.

"""
    _EDGES

Directed edges of the axes graph, as `(from, to)`. One editable home for the
connectivity; `_route` is derived from it, and `edge_theory` says which
theory's subgraph each belongs to.

Every edge here must have an `axes_rotation` method, and every such method must
appear here. The test suite checks both conditions.
"""
const _EDGES = (
    (ICRF(), MJ2000Eq()), (MJ2000Eq(), ICRF()),
    (MJ2000Eq(), MODEq()), (MODEq(), MJ2000Eq()),
    (MODEq(), TODEq()),    (TODEq(), MODEq()),
    (TODEq(), PEF()),      (PEF(), TODEq()),
    (PEF(), ITRF()),       (ITRF(), PEF()),
    (ICRF(), GCRF()),      (GCRF(), ICRF()),
    (GCRF(), CIRS()),      (CIRS(), GCRF()),
    (CIRS(), TIRS()),      (TIRS(), CIRS()),
    (TIRS(), ITRF()),      (ITRF(), TIRS()),
    (MJ2000Eq(), MJ2000Ec()), (MJ2000Ec(), MJ2000Eq()),
    (MODEq(), MODEc()),    (MODEc(), MODEq()),
    (TODEq(), TODEc()),    (TODEc(), TODEq()),
    (TODEq(), TEME()),     (TEME(), TODEq()),
)

# Frames that hang off a single hub are NOT listed here - they declare
# `hub_axes` instead, and their routes are derived from it. That covers the
# ones that cannot be listed (`CelestialBodyFixed{N}` is parametric; a user's
# frame does not exist yet) and, so there is one way to say it rather than two,
# the ones that could be.

# Breadth-first search over the edges available under `theory`, emitting a
# `_route` method for every pair that needs intermediates. Adjacent pairs are
# skipped: they already have a direct method.
# `nothing` as the theory means "use every edge", which generates the fallback
# set consulted when the active theory cannot reach the target on its own.
for theory in (FK5(), IAU2006(), nothing)
    allowed = theory === nothing ? _EDGES : filter(_EDGES) do (a, b)
        t = edge_theory(a, b)
        t === nothing || t === theory
    end

    nodes = unique(Iterators.flatten(allowed))

    for source in nodes
        # BFS from `source`, recording the predecessor of each node reached.
        prev    = Dict{Any,Any}()
        frontier = Any[source]
        while !isempty(frontier)
            next = Any[]
            for u in frontier, (a, b) in allowed
                a === u || continue
                (b === source || haskey(prev, b)) && continue
                prev[b] = u
                push!(next, b)
            end
            frontier = next
        end

        for target in nodes
            target === source && continue
            haskey(prev, target) || continue

            # Walk predecessors back to the source, then drop both endpoints:
            # `_route` carries only what lies between them.
            path = Any[target]
            while path[end] !== source
                push!(path, prev[path[end]])
            end
            between = Tuple(reverse(path)[2:end-1])
            isempty(between) && continue

            if theory === nothing
                @eval _route_any(::$(typeof(source)), ::$(typeof(target))) = $between
            else
                @eval _route(::$(typeof(source)), ::$(typeof(target)), ::$(typeof(theory))) =
                    $between
            end
        end

        # Reachability, including adjacent pairs. `_route` returns `nothing`
        # both for "adjacent, no intermediates needed" and for "no path at
        # all"; a derived route has to tell those apart, so record it here
        # where the search already knows.
        for target in nodes
            target === source && continue
            haskey(prev, target) || continue
            if theory === nothing
                @eval _reaches_any(::$(typeof(source)), ::$(typeof(target))) = true
            else
                @eval _reaches(::$(typeof(source)), ::$(typeof(target)), ::$(typeof(theory))) = true
            end
        end
    end
end

"""
    hub_axes(::Type{MyAxes}) -> axes or nothing

Returns the axes to which an extension frame connects.

Most extension frames are defined relative to one existing frame. Declaring
that frame allows AstroFrames to route transformations to all other reachable
frames:

```julia
struct PhobosFixed <: AstroFrames.AbstractAxes end
AstroFrames.hub_axes(::Type{PhobosFixed}) = ICRF()
```

# Arguments
- The axes type, not an instance; the result cannot depend on run-time data.

# Returns
The axes to which the extension frame is directly connected, or `nothing` for
a frame that belongs to a declared chain. The shipped Earth frames belong to
declared chains.

# Notes
An extension frame defines one `axes_rotation` method in each direction between
itself and its hub. No registry entry or transformation to other axes is
required. The declared hub must match the frame used by those two methods;
otherwise the extension frame remains unreachable.

# Example
```julia
struct PhobosFixed <: AstroFrames.AbstractAxes end

AstroFrames.hub_axes(::Type{PhobosFixed}) = ICRF()

function AstroFrames.axes_rotation(::ICRF, ::PhobosFixed, e::EpochScales)
    return body_axes_rotation(orientation_model(phobos), 401, epoch_tdb(e))
end
AstroFrames.axes_rotation(a::PhobosFixed, ::ICRF, e::EpochScales) =
    inv(axes_rotation(ICRF(), a, e))

axes_rotation(ITRF(), PhobosFixed(), epoch)    # routes, six edges, unassisted
```
"""
hub_axes(::Type{<:AbstractAxes}) = nothing

@inline hub_axes(a::AbstractAxes) = hub_axes(typeof(a))

# The frames that hang off ICRF. Declared here for the ones this package ships;
# a user declares their own the same way, which is the whole point.
hub_axes(::Type{<:CelestialBodyFixed}) = ICRF()
hub_axes(::Type{MoonPA})  = ICRF()
hub_axes(::Type{MoonME})  = ICRF()
hub_axes(::Type{RIC})     = ICRF()
hub_axes(::Type{LVLH})    = ICRF()
hub_axes(::Type{VNB})     = ICRF()

"""
    _reaches(source, target, theory) -> Bool

Whether a path exists, adjacent pairs included. Generated by the same search
that generates `_route`; `false` by default.
"""
_reaches(::AbstractAxes, ::AbstractAxes, ::AbstractFrameTheory) = false

"""
    _reaches_any(source, target) -> Bool

As `_reaches`, ignoring the active theory. The companion of `_route_any`.
"""
_reaches_any(::AbstractAxes, ::AbstractAxes) = false

# Adjacent-or-routed reachability for a pair the search enumerated.
@inline _connected(a, b, theory) = a === b || _reaches(a, b, theory)
@inline _connected_any(a, b) = a === b || _reaches_any(a, b)

"""
    _hub_route(source, target, theory) -> Tuple | Nothing

Route for a pair where at least one end hangs off a hub.

Derived rather than declared, so it covers frames the search could never have
enumerated — `CelestialBodyFixed{N}` is parametric over NAIF ID, and a user's
own frame does not exist when this package loads. `hub_axes` dispatches on the
type, so the whole thing resolves at compile time and nothing is searched on
the call path.
"""
@inline function _hub_route(source::AbstractAxes, target::AbstractAxes, theory)
    return _hub_route(source, target,
                      (a, b) -> _connected(a, b, theory),
                      (a, b) -> _route(a, b, theory))
end

"""
    _hub_route_any(source, target) -> Tuple | Nothing

As `_hub_route`, ignoring the active theory.

A frame connected through a hub inherits the hub's reachability. If reaching
the hub crosses the two Earth chains, reaching the extension frame does too,
and AstroFrames emits the `FR-FRAME-4` warning.
"""
@inline _hub_route_any(source::AbstractAxes, target::AbstractAxes) =
    _hub_route(source, target, _connected_any, _route_any)

@inline function _hub_route(source::AbstractAxes, target::AbstractAxes,
                            connected::F, route::G) where {F,G}
    hs = hub_axes(source)
    ht = hub_axes(target)

    if hs === nothing && ht === nothing
        return nothing                       # both interior; the search had its chance
    elseif hs !== nothing && ht !== nothing
        # Leaf to leaf. Same hub is one hop through it; different hubs go
        # through both, with whatever connects them in between.
        hs === ht && return (hs,)
        connected(hs, ht) || return nothing
        inner = route(hs, ht)
        return inner === nothing ? (hs, ht) : (hs, inner..., ht)
    elseif ht !== nothing
        # Interior source, leaf target: reach the hub, then step in.
        source === ht && return nothing      # adjacent; the direct edge handles it
        connected(source, ht) || return nothing
        inner = route(source, ht)
        return inner === nothing ? (ht,) : (inner..., ht)
    else
        # Leaf source, interior target: step out to the hub, then reach it.
        hs === target && return nothing      # adjacent
        connected(hs, target) || return nothing
        inner = route(hs, target)
        return inner === nothing ? (hs,) : (hs, inner...)
    end
end

# Recursive composition over the declared tuple. Written this way rather than
# as a loop so the tuple length is known to the compiler and the whole chain
# inlines into a fixed sequence of matrix products.
@inline _compose(from::AbstractAxes, to::AbstractAxes, ::Tuple{}, e::EpochScales) =
    axes_rotation(from, to, e)

@inline _compose(from::AbstractAxes, to::AbstractAxes, r::Tuple, e::EpochScales) =
    _compose(first(r), to, Base.tail(r), e) * axes_rotation(from, first(r), e)

# `frame_theory()` reads a `Ref{AbstractFrameTheory}`, so its return type is
# abstract. The barrier confines that to a single call: `theory` is concrete
# inside `_routed_rotation`, so `_route` resolves statically and `_compose`
# inlines.
#
# Measured, this is not where the time goes — a routed FK5 chain costs 3.81 µs
# against 3.74 µs hand-composed, so routing is free either way. The barrier is
# here because an abstract global should not reach a dispatch site, not
# because it fixed a measured problem.
function axes_rotation(source::A1, target::A2, e::EpochScales) where {A1<:AbstractAxes, A2<:AbstractAxes}
    return _routed_rotation(source, target, e, frame_theory())
end

function _routed_rotation(source::A1, target::A2, e::EpochScales,
                          theory::T) where {A1<:AbstractAxes, A2<:AbstractAxes, T<:AbstractFrameTheory}
    route = _route(source, target, theory)
    route === nothing || return _compose(source, target, route, e)

    # Frames that hang off a hub are not in the search's node list — a user's
    # own frame did not exist when it ran — so their routes are derived.
    hub = _hub_route(source, target, theory)
    hub === nothing || return _compose(source, target, hub, e)

    # No route within the active theory. A route may still exist across the
    # two chains; that is a modelling choice rather than an error
    # (`FR-FRAME-4`), so it proceeds and warns once.
    fallback = _route_any(source, target)
    if fallback !== nothing
        _warn_out_of_theory(source, target, theory)
        return _compose(source, target, fallback, e)
    end

    hub_fallback = _hub_route_any(source, target)
    if hub_fallback !== nothing
        _warn_out_of_theory(source, target, theory)
        return _compose(source, target, hub_fallback, e)
    end

    throw(ArgumentError(
        "No state transform registered from $(A1) to $(A2). " *
        "Supported edges in this release: identity, ICRF ↔ CelestialBodyFixed{N}, " *
        "ICRF ↔ MoonPA, ICRF ↔ MoonME, ICRF ↔ MJ2000Eq, MJ2000Eq ↔ MJ2000Ec, " *
        "MJ2000Eq ↔ MODEq, MODEq ↔ TODEq, TODEq ↔ PEF, PEF ↔ ITRF, " *
        "MODEq ↔ MODEc, TODEq ↔ TODEc, TODEq ↔ TEME, " *
        "ICRF ↔ GCRF, GCRF ↔ CIRS, CIRS ↔ TIRS, TIRS ↔ ITRF. " *
        "ICRF ↔ RIC, ICRF ↔ LVLH, ICRF ↔ VNB. " *
        "A frame of your own needs an `axes_rotation` method and a `hub_axes` " *
        "declaration saying which frame it attaches to; see `hub_axes`."))
end

# --- The same routing, carrying parameters -----------------------------------
#
# Some frames are defined by data this package cannot obtain by reaching
# downward — an orbit-relative frame is defined by another object's
# trajectory, which lives above AstroFrames. `params` is how the caller
# supplies it, and it has to survive composition: `ITRF → VNB` routes through
# five edges, and only the last one consumes it.
#
# The rule is that an edge needing no parameters ignores them. A caller can
# therefore pass `params` uniformly without knowing which edges along a route
# will use it, and an edge that does need them wins by dispatch specificity.

@inline _compose(from::AbstractAxes, to::AbstractAxes, ::Tuple{},
                 e::EpochScales, p::NamedTuple) = axes_rotation(from, to, e, p)

@inline _compose(from::AbstractAxes, to::AbstractAxes, r::Tuple,
                 e::EpochScales, p::NamedTuple) =
    _compose(first(r), to, Base.tail(r), e, p) * axes_rotation(from, first(r), e, p)

# The parameter-carrying form documented with `axes_rotation` above.
function axes_rotation(source::A1, target::A2, e::EpochScales,
                       p::NamedTuple) where {A1<:AbstractAxes, A2<:AbstractAxes}
    return _routed_rotation(source, target, e, p, frame_theory())
end

function _routed_rotation(source::A1, target::A2, e::EpochScales, p::NamedTuple,
                          theory::T) where {A1<:AbstractAxes, A2<:AbstractAxes, T<:AbstractFrameTheory}
    route = _route(source, target, theory)
    route === nothing || return _compose(source, target, route, e, p)

    hub = _hub_route(source, target, theory)
    hub === nothing || return _compose(source, target, hub, e, p)

    fallback = _route_any(source, target)
    if fallback !== nothing
        _warn_out_of_theory(source, target, theory)
        return _compose(source, target, fallback, e, p)
    end

    hub_fallback = _hub_route_any(source, target)
    if hub_fallback !== nothing
        _warn_out_of_theory(source, target, theory)
        return _compose(source, target, hub_fallback, e, p)
    end

    # Adjacent, or genuinely unsupported. Either way the parameterless form is
    # the right thing to call: it holds the direct edges, and it raises the
    # error naming what is missing when there is no edge at all.
    return axes_rotation(source, target, e)
end

# --- Helper: invert a state-transform matrix -------------------------------
#
# For a state-transform of the SPICE block form `[R 0; Ṙ R]` where R is a
# rotation matrix (orthogonal), the inverse is `[Rᵀ 0; Ṙᵀ Rᵀ]`. This is
# faster and more numerically stable than a generic 6×6 inversion.

function _invert_rotation(M::AbstractMatrix)
    R  = @SMatrix [M[1,1] M[1,2] M[1,3]
                   M[2,1] M[2,2] M[2,3]
                   M[3,1] M[3,2] M[3,3]]
    Ṙ  = @SMatrix [M[4,1] M[4,2] M[4,3]
                   M[5,1] M[5,2] M[5,3]
                   M[6,1] M[6,2] M[6,3]]
    Rᵀ = transpose(R)
    Ṙᵀ = transpose(Ṙ)
    Z  = @SMatrix zeros(3, 3)
    return @SMatrix [Rᵀ[1,1] Rᵀ[1,2] Rᵀ[1,3] Z[1,1] Z[1,2] Z[1,3]
                     Rᵀ[2,1] Rᵀ[2,2] Rᵀ[2,3] Z[2,1] Z[2,2] Z[2,3]
                     Rᵀ[3,1] Rᵀ[3,2] Rᵀ[3,3] Z[3,1] Z[3,2] Z[3,3]
                     Ṙᵀ[1,1] Ṙᵀ[1,2] Ṙᵀ[1,3] Rᵀ[1,1] Rᵀ[1,2] Rᵀ[1,3]
                     Ṙᵀ[2,1] Ṙᵀ[2,2] Ṙᵀ[2,3] Rᵀ[2,1] Rᵀ[2,2] Rᵀ[2,3]
                     Ṙᵀ[3,1] Ṙᵀ[3,2] Ṙᵀ[3,3] Rᵀ[3,1] Rᵀ[3,2] Rᵀ[3,3]]
end

# --- Public entry points -----------------------------------------------------
#
# One conversion, at the boundary. Everything above this line works in plain
# numbers.

axes_rotation(from::AbstractAxes, to::AbstractAxes, t::Time) =
    axes_rotation(from, to, _scales(t))

axes_rotation(from::AbstractAxes, to::AbstractAxes, jd_tdb::Real) =
    axes_rotation(from, to, _scales(jd_tdb))

axes_rotation(from::AbstractAxes, to::AbstractAxes, t::Time, p::NamedTuple) =
    axes_rotation(from, to, _scales(t), p)

axes_rotation(from::AbstractAxes, to::AbstractAxes, jd_tdb::Real, p::NamedTuple) =
    axes_rotation(from, to, _scales(jd_tdb), p)
