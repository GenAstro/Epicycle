# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Orbit-relative axes.
#
# The family whose defining data cannot be obtained downward. Earth frames read
# EOP and body-fixed frames read IAU constants, both of which AstroUniverse
# holds. An orbit-relative frame is defined by *another object's trajectory*,
# which lives above AstroFrames — so the caller passes it in:
#
#     axes_rotation(ICRF(), RIC(), epoch, (; reference_state))
#     axes_rotation(ICRF(), VNB(), epoch, (; reference_state, reference_accel))
#
# `params` carries evaluated values, never a spacecraft or a force model.
# AstroFrames never reaches up; the caller reaches down.
#
# All three are one construction — a primary direction, a secondary direction,
# and a labelling — differing only in those choices:
#
#     RIC   primary r̂,  secondary ĥ    radial, in-track, cross-track
#     LVLH  primary −r̂, secondary −ĥ   nadir-pointing
#     VNB   primary v̂,  secondary ĥ    velocity, normal, binormal
#
# **Radial-primary frames need only a state.** `r̂` turns at `ω = h/r²`, which
# is a function of position and velocity alone. Velocity-primary frames do not:
# `v̂` turns at a rate set by the force environment. With no acceleration,
# VNB describes an unaccelerated reference and has zero rate. RIC and LVLH
# still turn with the radius, while omitting orbit-plane acceleration.
# =============================================================================

function _orbit_relative_missing(frame, field, why)
    return ArgumentError(
        "$(frame) axes need `$(field)`: $(why). Make a spacecraft the coordinate system's " *
        "origin, `CoordinateSystem(chief, $(frame)())`, and its state is used; or pass the " *
        "reference orbit yourself, `Coordinate(c, cs, (; reference_state = x))` with `x` the " *
        "reference position and velocity in ICRF. AstroFrames does not propagate.")
end

# A key that starts `reference` but is not one of the two read here is almost certainly a
# misspelling, and ignoring it silently left a frame without the acceleration its author meant to
# give. Other keys pass: params travel the whole route, and a frame of the user's own may read
# keys of its own.
function _check_reference_keys(p::NamedTuple, frame)
    for k in keys(p)
        k in (:reference_state, :reference_accel) && continue
        startswith(String(k), "reference") && throw(ArgumentError(
            "$(frame) axes: unknown parameter `$(k)`. The reference orbit is given by " *
            "`reference_state` and, optionally, `reference_accel`."))
    end
    return nothing
end

# A reference orbit at the origin, or moving along its own radius, has no plane, and the frame's
# axes come out NaN. Kept as NaN, as the state conversions do, so an optimiser can step back; but
# said once, since otherwise nothing names the cause.
@noinline function _warn_degenerate_reference(frame)
    @warn "AstroFrames: the reference orbit for $(frame) axes has zero angular momentum (zero " *
          "position, zero velocity, or velocity along the radius), so the axes are undefined " *
          "and come out NaN. Shown once per session." maxlog = 1
    return nothing
end

"""
    _unit_and_rate(u, u̇) -> (û, d(û)/dt)

Unit vector and its time derivative.

Differentiating `u/|u|` gives `(u̇ − û(û·u̇))/|u|` — the component of `u̇`
perpendicular to `u`, scaled by `1/|u|`. The parallel component changes the
vector's length, not its direction, so it drops out.
"""
@inline function _unit_and_rate(u::AbstractVector, u̇::AbstractVector)
    n = norm(u)
    û = u / n
    return û, (u̇ - û * dot(û, u̇)) / n
end

"""
    _triad(primary, primary_rate, secondary, secondary_rate) -> (R, Ṙ)

Build an orthonormal triad and its rate from a primary and secondary direction.

Rows of `R` are the target-frame axes expressed in the source frame, in the
order (primary, secondary, completion), so `R * r_source` gives the components
along them. The third axis is `primary × secondary`, which makes the set
right-handed in that order.
"""
@inline function _triad(p::AbstractVector, ṗ::AbstractVector,
                        s::AbstractVector, ṡ::AbstractVector)
    c  = cross(p, s)
    ċ  = cross(ṗ, s) + cross(p, ṡ)
    R  = @SMatrix [p[1] p[2] p[3]; s[1] s[2] s[3]; c[1] c[2] c[3]]
    Ṙ  = @SMatrix [ṗ[1] ṗ[2] ṗ[3]; ṡ[1] ṡ[2] ṡ[3]; ċ[1] ċ[2] ċ[3]]
    return R, Ṙ
end

"""
    needs_reference_orbit(axes) -> Bool

Whether these axes are defined by a reference orbit rather than by a body.

# Arguments
- `axes::AbstractAxes` — the axes to ask about, such as `RIC()`.

# Returns
`true` for [`RIC`](@ref), [`LVLH`](@ref) and [`VNB`](@ref), which need a
reference state to define their orientation, and `false` for axes fixed to a
body or to the sky.

# Notes
The default is `false`. A user-defined orbit-relative frame overrides this
method to receive automatic reference-orbit handling in coordinate conversions.

# Example

```jldoctest
needs_reference_orbit(RIC()), needs_reference_orbit(ICRF())

# output
(true, false)
```
"""
needs_reference_orbit(::AbstractAxes) = false
needs_reference_orbit(::RIC)  = true
needs_reference_orbit(::LVLH) = true
needs_reference_orbit(::VNB)  = true

"""
    _reference(p::NamedTuple, frame) -> (r, v, a)

Pull the reference orbit out of `params`, failing by name if it is absent.

`a` defaults to zero when acceleration is not supplied. For a radial-primary
frame that is exact under two-body motion and neglects only the rotation of
the orbit plane under out-of-plane forces. Without supplied acceleration, VNB describes an unaccelerated reference with zero rate.
"""
@inline function _reference(p::NamedTuple, frame)
    _check_reference_keys(p, frame)
    haskey(p, :reference_state) || throw(_orbit_relative_missing(
        frame, :reference_state, "the frame is defined by a reference orbit"))
    s = p.reference_state
    s isa Union{AbstractVector, Tuple} || throw(ArgumentError(
        "$(frame) axes: `reference_state` must contain six real numbers; got $(repr(s))."))
    length(s) == 6 || throw(ArgumentError(
        "$(frame) axes: `reference_state` is a position and velocity, six numbers in km and " *
        "km/s; got $(length(s))."))
    all(x -> x isa Real, s) || throw(ArgumentError(
        "$(frame) axes: `reference_state` must contain six real numbers; got $(s)."))
    r = SVector{3}(s[1], s[2], s[3])
    v = SVector{3}(s[4], s[5], s[6])
    if haskey(p, :reference_accel)
        acc = p.reference_accel
        acc isa Union{AbstractVector, Tuple} || throw(ArgumentError(
            "$(frame) axes: `reference_accel` must contain three real numbers; got $(repr(acc))."))
        length(acc) == 3 || throw(ArgumentError(
            "$(frame) axes: `reference_accel` is an acceleration, three numbers in km/s²; got " *
            "$(length(acc))."))
        all(x -> x isa Real, acc) || throw(ArgumentError(
            "$(frame) axes: `reference_accel` must contain three real numbers; got $(acc)."))
        a = SVector{3}(acc[1], acc[2], acc[3])
    else
        a = zero(SVector{3,Float64})
    end
    # Squaring before the zero check also works for dual inputs: at zero angular momentum the
    # squared magnitude and its first derivatives are zero, even when the vector's are not.
    h2 = sum(abs2, cross(r, v))
    if any(isnan, s) || any(isnan, a)
        bad = (r + v + a) .* NaN
        return bad, bad, bad
    elseif !all(isfinite, s) || !all(isfinite, a) || !isfinite(h2) ||
           !all(isfinite, cross(r, a))
        @warn "AstroFrames: the reference data for $(frame) axes are non-finite or overflow; the axes are undefined and return NaNs." maxlog = 1
        bad = (r + v + a) .* NaN
        return bad, bad, bad
    elseif !(h2 > 0)
        _warn_degenerate_reference(frame)
        bad = (r + v + a) .* NaN
        return bad, bad, bad
    end
    return r, v, a
end

# --- RIC: radial, in-track, cross-track -------------------------------------

"""
    axes_rotation(::ICRF, ::RIC, epoch, params) -> SMatrix{6,6}

Radial / in-track / cross-track axes of a reference orbit, in that order.

`R̂` follows the radius, `Ĉ` follows the orbit normal `r × v`, and `Î = Ĉ × R̂`
completes the set. The `Î` axis is transverse to the radius in the orbit plane;
it need not align with velocity on an eccentric orbit. The rows are `(R̂, Î, Ĉ)`:
radial, in-track, and cross-track components.

`params` must carry `reference_state`; `reference_accel` is optional and, when
given, includes the rotation of the orbit plane under out-of-plane forces.

Also called RSW or QSW; AstroFrames provides this convention as `RIC`.
"""
function axes_rotation(::ICRF, ::RIC, e::EpochScales, p::NamedTuple)
    r, v, a = _reference(p, "RIC")
    R̂, R̂̇ = _unit_and_rate(r, v)
    Ĉ, Ĉ̇ = _unit_and_rate(cross(r, v), cross(r, a))
    # `_triad` returns rows (primary, secondary, primary × secondary). Passing
    # (Ĉ, R̂) makes the third row Ĉ × R̂ = Î, the in-track direction; taking
    # (R̂, Ĉ) instead would put −Î there, pointing against the motion.
    R, Ṙ = _triad(Ĉ, Ĉ̇, R̂, R̂̇)
    # Rows arrive as (Ĉ, R̂, Î); RIC wants (R̂, Î, Ĉ).
    P = @SMatrix [0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0]
    return _rotation_with_rate(P * R, P * Ṙ)
end

# --- LVLH: local vertical, local horizontal ---------------------------------

"""
    axes_rotation(::ICRF, ::LVLH, epoch, params) -> SMatrix{6,6}

Local-vertical / local-horizontal axes of a reference orbit.

Uses the nadir-pointing convention: `ẑ = −r̂`, `ŷ = −ĥ` is opposite the orbit
normal, and `x̂ = ŷ × ẑ` completes the set and lies nearly along-track.
LVLH definitions vary in axis order and signs; check the convention when
exchanging data with another system.

`params` must carry `reference_state`; `reference_accel` is optional.
"""
function axes_rotation(::ICRF, ::LVLH, e::EpochScales, p::NamedTuple)
    r, v, a = _reference(p, "LVLH")
    r̂, r̂̇ = _unit_and_rate(r, v)
    ĥ, ĥ̇ = _unit_and_rate(cross(r, v), cross(r, a))
    # ŷ = −ĥ, ẑ = −r̂, x̂ = ŷ × ẑ. Ordering the triad as (x̂, ŷ, ẑ) means
    # passing ŷ and ẑ as primary and secondary, since ŷ × ẑ = x̂.
    R, Ṙ = _triad(-ĥ, -ĥ̇, -r̂, -r̂̇)
    # `_triad` returns rows (primary, secondary, completion) = (ŷ, ẑ, x̂);
    # LVLH wants (x̂, ŷ, ẑ), so rotate the rows.
    P = @SMatrix [0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0]
    return _rotation_with_rate(P * R, P * Ṙ)
end

# --- VNB: velocity, normal, binormal ----------------------------------------

"""
    axes_rotation(::ICRF, ::VNB, epoch, params) -> SMatrix{6,6}

Velocity / normal / binormal axes of a reference orbit.

`V̂` along velocity, `N̂` along the orbit normal, `B̂ = V̂ × N̂`.

`params` must carry `reference_state` in ICRF, in km and km/s. The optional
`reference_accel`, in km/s², determines the rate of the velocity direction and
orbit normal. Without it, the frame rate is zero, representing an
unaccelerated reference. This is suitable for resolving instantaneous vector
components. Supply acceleration to include the velocity contribution from
the frame's motion.
"""
function axes_rotation(::ICRF, ::VNB, e::EpochScales, p::NamedTuple)
    r, v, a = _reference(p, "VNB")
    V̂, V̂̇ = _unit_and_rate(v, a)
    N̂, N̂̇ = _unit_and_rate(cross(r, v), cross(r, a))
    return _rotation_with_rate(_triad(V̂, V̂̇, N̂, N̂̇)...)
end

# --- Reverse directions ------------------------------------------------------

for F in (:RIC, :LVLH, :VNB)
    @eval axes_rotation(::$F, ::ICRF, e::EpochScales, p::NamedTuple) =
        _invert_rotation(axes_rotation(ICRF(), $F(), e, p))
end

# --- Calling without parameters at all ---------------------------------------
#
# Must say what is missing rather than failing as a `MethodError` several calls
# down.

for F in (:RIC, :LVLH, :VNB)
    @eval begin
        axes_rotation(::ICRF, ::$F, ::EpochScales) = throw(_orbit_relative_missing(
            $(string(F)), :reference_state, "the frame is defined by a reference orbit"))
        axes_rotation(::$F, ::ICRF, ::EpochScales) = throw(_orbit_relative_missing(
            $(string(F)), :reference_state, "the frame is defined by a reference orbit"))
    end
end
