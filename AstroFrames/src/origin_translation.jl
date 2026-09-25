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

Offset to add after rotating, when the origin changes.

Returns the **state of `source_origin` relative to `target_origin`, expressed
in `axes`** — position in km, velocity in km/s. Both halves of that sentence
are load-bearing:

  * **Direction.** It is the source origin seen from the target, not the
    reverse. Check it against the degenerate case: a point sitting at the
    source origin has zero coordinates there, so its coordinates about the
    target origin must be exactly this offset.
  * **Axes.** An offset is a vector, and a vector expressed in the wrong axes
    is wrong by a rotation while looking entirely plausible. Nothing about the
    numbers reveals the mistake, so the axes are named rather than assumed.

Composed with a rotation, the full frame change is

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

# Notes
- The ephemeris is delivered in ICRF and is rotated into `axes` here, so a
  caller working in an FK5 frame gets the frame bias applied automatically
  rather than having to remember it. Forgetting it is ≈16 km at 1 AU.
- **Not differentiable**: the ephemeris comes from SPICE. An origin shift is
  therefore outside any AD path.

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
    own = axes_rotation(cs.axes, ICRF(), e) * SVector{6}(_cartesian_vector(origin))
    body, rest = _resolve_origin(cs.origin, e)
    return (body, own + rest)
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
        between = SVector{6}(translate_state(target_body, source_body, e.tdb))
        return axes_rotation(ICRF(), axes, e, params) *
               (source_offset + between - target_offset)
    end

    # `translate_state(from, to)` is the state of `to` seen from `from`, so the
    # arguments are reversed here: we want the SOURCE origin seen from the
    # TARGET origin.
    Δ_icrf = translate_state(target_origin, source_origin, e.tdb)

    # The ephemeris is ICRF; rotate it into the requested axes. Tagging it ICRF
    # rather than passing SPICE's own "J2000" label through is what makes the
    # frame bias an ordinary graph edge instead of something to remember.
    #
    # `params` matters when the target axes are defined by a reference orbit:
    # the offset has to land in the same axes the state was rotated into.
    return axes_rotation(ICRF(), axes, e, params) * SVector{6}(Δ_icrf)
end

# --- Public entry points -----------------------------------------------------

origin_translation(from::AbstractPoint, to::AbstractPoint, axes::AbstractAxes, t::Time,
                   params::NamedTuple = NamedTuple()) =
    origin_translation(from, to, axes, _scales(t), params)

origin_translation(from::AbstractPoint, to::AbstractPoint, axes::AbstractAxes, jd_tdb::Real,
                   params::NamedTuple = NamedTuple()) =
    origin_translation(from, to, axes, _scales(jd_tdb), params)
