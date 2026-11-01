# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Origin translation — the second primitive.
#
# A frame change is two *kinds* of operation, and this is the other one. A
# rotation is linear; a translation is not, so forcing both into one object
# would mean a 7×7 homogeneous form that buys nothing.
#
#     state_target = axes_rotation(...) * state_source + origin_translation(...)
#
# In general the two alternate rather than happening once each: Earth-fixed to
# a rotating libration-point frame is rotate, translate, rotate. The
# translation happens in ICRF, because that is the frame the ephemeris arrives
# in — see the `axes` argument below.
# =============================================================================

using AstroUniverse: translate_state

"""
    origin_translation(source_origin, target_origin, axes, epoch[, params]) -> SVector{6}

Return the state of `source_origin` relative to `target_origin`, expressed in
`axes`: position in km followed by velocity in km/s.

Add this offset after rotating a state into the target axes. A point at the
source origin has zero coordinates there, so its target coordinates equal
this offset. For a frame change without additional reference parameters:

    state_target = axes_rotation(source_axes, target_axes, epoch) * state_source
                 + origin_translation(source_origin, target_origin, target_axes, epoch)

# Arguments
- `source_origin::AbstractPoint` — the origin the state is currently about, such
  as `earth`.
- `target_origin::AbstractPoint` — the origin to move it to, such as `moon`.
- `axes::AbstractAxes` — the axes the offset is expressed in; pass the target axes.
- `epoch` — an `AstroEpochs.Time`, which carries its own time scale, or a Julian
  date in TDB.
- `params::NamedTuple` — optional; the same evaluated numbers
  [`axes_rotation`](@ref) takes, needed when `axes` are defined by a reference
  orbit.

# Returns
An `SVector{6}` containing the source origin's position and velocity relative
to the target, in the requested axes.

# Notes
Ephemeris states arrive in ICRF and are rotated into `axes`, including the
frame bias when required. The bias corresponds to about 16 km at 1 AU.

SPICE ephemeris calls do not support epoch differentiation. The returned
offset can still be added to a differentiable state at a fixed ordinary epoch.

# Example

```jldoctest
using AstroUniverse: earth, moon
using AstroEpochs: Time, UTC, ISOT
Δ = origin_translation(earth, moon, ICRF(), Time("2020-01-01T00:00:00.000", UTC(), ISOT()))
length(Δ)

# output
6
```
"""
function origin_translation end

"""
    _resolve_origin(origin, e) -> (body, offset)

An origin as a body the ephemeris knows, plus an ICRF offset from it.

A body resolves to itself with no offset. An origin that carries its own state,
such as a spacecraft, resolves to the body at the center of its frame plus that
state. The recursion ends at a body because a carried state is expressed about
one.

Returning the body rather than a state relative to some reference is what lets
*both* origins carry their own state: each resolves independently, and the
ephemeris is consulted once, between the two bodies.
"""
function _resolve_origin(origin, e::EpochScales)
    carries_own_state(origin) || return (origin, zero(SVector{6,Float64}))

    _require_epoch_match(origin, e.tdb,
        "$(nameof(typeof(origin))) is the origin and carries its own state")

    cs  = frame_of(origin)
    params = _reference_orbit_from_origin(cs, e, NamedTuple())
    own = axes_rotation(cs.axes, ICRF(), e, params) * SVector{6}(_cartesian_vector(origin))
    body, rest = _resolve_origin(cs.origin, e)
    return (body, own + rest)
end

# The TDB date for the ephemeris. SPICE takes a Float64, so a dual date, carrying a derivative,
# cannot pass; say so, as the lunar axes do, rather than fail inside the conversion.
@inline function _ephemeris_date(e::EpochScales)
    t = e.tdb
    t isa Union{AbstractFloat, Integer} || throw(ArgumentError(
        "A change of origin reads the ephemeris from SPICE, which is not differentiable, so it " *
        "cannot take an epoch of type $(typeof(t)). Differentiate with respect to the state " *
        "instead, or keep the origin unchanged."))
    return t
end

function origin_translation(source_origin::AbstractPoint,
                            target_origin::AbstractPoint,
                            axes::AbstractAxes,
                            e::EpochScales,
                            params::NamedTuple = NamedTuple())
    # An origin that carries its own state is not in the ephemeris, so the
    # offset is built from what it carries instead. Kept off the body-to-body
    # path, which is the common one.
    if carries_own_state(source_origin) || carries_own_state(target_origin)
        source_body, source_offset = _resolve_origin(source_origin, e)
        target_body, target_offset = _resolve_origin(target_origin, e)

        # Both measured from the target's body, then differenced: the source
        # origin seen from the target origin.
        between = SVector{6}(translate_state(target_body, source_body, _ephemeris_date(e)))
        return axes_rotation(ICRF(), axes, e, params) *
               (source_offset + between - target_offset)
    end

    # `translate_state(from, to)` is the state of `to` seen from `from`, so the
    # arguments are reversed here: we want the SOURCE origin seen from the
    # TARGET origin.
    Δ_icrf = translate_state(target_origin, source_origin, _ephemeris_date(e))

    # The ephemeris is ICRF; rotate it into the requested axes. Tagging it ICRF
    # rather than passing SPICE's own "J2000" label through is what makes the
    # frame bias an ordinary graph edge instead of something to remember.
    #
    # `params` matters when the target axes are defined by a reference orbit:
    # the offset has to land in the same axes the state was rotated into.
    return axes_rotation(ICRF(), axes, e, params) * SVector{6}(Δ_icrf)
end

# --- Public entry points -----------------------------------------------------

function origin_translation(from::AbstractPoint, to::AbstractPoint, axes::AbstractAxes,
                            epoch::Union{Time, Real}, params::NamedTuple = NamedTuple())
    _valid_epoch(epoch) || return SVector{6}(_nan_epoch_values(epoch, Val(6)))
    return origin_translation(from, to, axes, _scales(epoch), params)
end
origin_translation(from::AbstractPoint, to::AbstractPoint, axes::AbstractAxes, epoch) =
    throw(ArgumentError("origin_translation: epoch must be an AstroEpochs.Time or a TDB Julian date; got $(repr(epoch)) ($(typeof(epoch)))."))

origin_translation(from::AbstractPoint, to::AbstractPoint, axes::AbstractAxes, epoch, params) =
    throw(ArgumentError("origin_translation: params must be a NamedTuple of evaluated reference data; got $(repr(params)) ($(typeof(params)))."))
origin_translation(from::AbstractPoint, to::AbstractPoint, axes::AbstractAxes, epoch, params::NamedTuple) =
    origin_translation(from, to, axes, epoch)
