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
# `v̂` turns at a rate set by the force environment, so `VNB` needs
# acceleration and says so. This is the whole practical difference between
# them, and it is why `reference_accel` is required for one and optional for
# the others.
# =============================================================================

function _orbit_relative_missing(frame, field, why)
    return ArgumentError(
        "$(frame) axes need `$(field)` in params: $(why). " *
        "Pass it as `axes_rotation(source, $(frame)(), epoch, (; reference_state))`. " *
        "The reference orbit comes from the caller; AstroFrames does not propagate.")
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
reference state before they mean anything, and `false` for axes fixed to a body
or to the sky.

# Example

```jldoctest
needs_reference_orbit(RIC()), needs_reference_orbit(ICRF())

# output
(true, false)
```

`false` by default. A user frame built from a reference orbit says so by
adding a method, and then gets the same automatic handling the shipped
orbit-relative frames do.
"""
needs_reference_orbit(::AbstractAxes) = false
needs_reference_orbit(::RIC)  = true
needs_reference_orbit(::LVLH) = true
needs_reference_orbit(::VNB)  = true

"""
    _reference(p::NamedTuple, frame, needs_accel::Bool) -> (r, v, a)

Pull the reference orbit out of `params`, failing by name if it is absent.

`a` defaults to zero when acceleration is not supplied. For a radial-primary
frame that is exact under two-body motion and neglects only the rotation of
the orbit plane under out-of-plane forces. Velocity-primary frames require
acceleration to determine their rate.
"""
@inline function _reference(p::NamedTuple, frame)
    haskey(p, :reference_state) || throw(_orbit_relative_missing(
        frame, :reference_state, "the frame is defined by a reference orbit"))
    s = p.reference_state
    r = SVector{3}(s[1], s[2], s[3])
    v = SVector{3}(s[4], s[5], s[6])
    a = haskey(p, :reference_accel) ? SVector{3}(p.reference_accel) : zero(SVector{3,Float64})
    return r, v, a
end

# --- RIC: radial, in-track, cross-track -------------------------------------

"""
    axes_rotation(::ICRF, ::RIC, epoch, params) -> SMatrix{6,6,Float64,36}

Radial / in-track / cross-track axes of a reference orbit, in that order.

`R̂` follows the radius, `Ĉ` follows the orbit normal `r × v`, and `Î = Ĉ × R̂`
completes the set. The `Î` axis is exactly along-track only for a circular orbit. The
rows come back as `(R̂, Î, Ĉ)`, so the second component of a transformed vector
is in-track and the third is cross-track, which is what the name says and what
every other tool means by it.

`params` must carry `reference_state`; `reference_accel` is optional and, when
given, includes the rotation of the orbit plane under out-of-plane forces.

Identical to the frames named RSW and QSW elsewhere; those are not separate
frames and are not provided under separate names.
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
    axes_rotation(::ICRF, ::LVLH, epoch, params) -> SMatrix{6,6,Float64,36}

Local-vertical / local-horizontal axes of a reference orbit.

**Convention, stated because it varies by source:** `ẑ` is nadir, `−r̂`; `ŷ` is
the negative orbit normal, `−ĥ`; `x̂ = ŷ × ẑ` completes the set and lies nearly
along-track. This is the nadir-pointing convention used for spacecraft
attitude. Sources differ on the signs, so anything comparing against another
tool should check theirs rather than assume.

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
    axes_rotation(::ICRF, ::VNB, epoch, params) -> SMatrix{6,6,Float64,36}

Velocity / normal / binormal axes of a reference orbit.

`V̂` along velocity, `N̂` along the orbit normal, `B̂ = V̂ × N̂`.

`params` must carry `reference_state`. `reference_accel` is optional, and what
it changes is worth understanding, because VNB turns for a different reason
than [`RIC`](@ref) and [`LVLH`](@ref) do.

Those frames are built from the radius. The radius direction turns with the
orbit under any force model, and its rate follows from the supplied velocity.

VNB is built from the velocity. A velocity direction turns only when the
velocity changes, and velocity changes only under acceleration. So VNB's rate
is not a property of the geometry the way RIC's is; it is a property of the
force environment.

Acceleration is therefore a separate input. Without `reference_accel`, the
velocity direction and orbit normal remain fixed and the frame rate is zero.
This describes an unaccelerated reference and supports instantaneous resolution
into along-track, normal, and binormal components. Supplying acceleration
includes the frame rate of an accelerating reference.
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
