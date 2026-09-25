# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Rotation-matrix builders.
#
# What the Earth chain needs, and the three ways a 3×3 rotation becomes the
# 6×6 that carries velocity with it:
#
#         ⎡ R    0 ⎤
#   M  =  ⎢        ⎥       state_target = M · state_source
#         ⎣ Ṙ    R ⎦
#
# `Ṙ` is what makes the difference between rotating a position and rotating a
# state — it carries the ω × r term. A frame with no rate uses a zero block, a
# frame spinning about one axis gives its rate as a scalar, and anything else
# supplies Ṙ directly.
#
# The pole-and-prime-meridian construction that used to live here is now
# `AstroUniverse.pole_axes_rotation`, because how a body is oriented belongs
# to the body.
# =============================================================================

# --- Elemental passive rotations --------------------------------------------
#
# X for the obliquity edges, Z for Earth rotation in both theories.

@inline function _Rz(θ::Real)
    c, s = cos(θ), sin(θ)
    return @SMatrix [ c  s  0.0
                     -s  c  0.0
                     0.0 0.0 1.0]
end

@inline function _Rx(θ::Real)
    c, s = cos(θ), sin(θ)
    return @SMatrix [1.0 0.0 0.0
                     0.0  c   s
                     0.0 -s   c]
end

# The pole-and-prime-meridian construction that used to live here now lives in
# `AstroUniverse` as `pole_axes_rotation`, because how a body is oriented is a
# property of the body and models are registered there. This file keeps the
# rotations the Earth chain needs.

"""
    _rotation_no_rate(R::AbstractMatrix) -> SMatrix{6,6,Float64,36}

Assemble the 6×6 state rotation for a 3×3 rotation `R` whose rate is zero:

        ⎡ R   0 ⎤
    M = ⎢       ⎥
        ⎣ 0   R ⎦

Velocity rotates by `R` with no rate coupling. Two kinds of edge use this:

  * genuinely static ones, including frame bias and the J2000 obliquity, where
    `Ṙ` is exactly zero;
  * precession, nutation, and polar-motion edges, where `Ṙ` is *neglected*
    per the rate policy: their true rates are ~1e-11 rad/s against Earth
    spin at 7.3e-5, and this is what SatelliteToolboxTransformations'
    own velocity handling assumes.

**Never use it for a spin edge** (Earth rotation via ERA or GAST). There the
rate is the dominant term, and omitting it corrupts velocity while leaving
position unchanged.
"""
function _rotation_no_rate(R::AbstractMatrix)
    Z = @SMatrix zeros(3, 3)
    return @SMatrix [R[1,1] R[1,2] R[1,3] Z[1,1] Z[1,2] Z[1,3]
                     R[2,1] R[2,2] R[2,3] Z[2,1] Z[2,2] Z[2,3]
                     R[3,1] R[3,2] R[3,3] Z[3,1] Z[3,2] Z[3,3]
                     Z[1,1] Z[1,2] Z[1,3] R[1,1] R[1,2] R[1,3]
                     Z[2,1] Z[2,2] Z[2,3] R[2,1] R[2,2] R[2,3]
                     Z[3,1] Z[3,2] Z[3,3] R[3,1] R[3,2] R[3,3]]
end

"""
    _rotation_with_spin(R::AbstractMatrix, ω::Real) -> SMatrix{6,6,Float64,36}

Assemble the 6×6 state rotation for an edge whose target frame **rotates**
about its own z-axis at rate `ω` [rad/s] relative to the source:

        ⎡  R       0 ⎤
    M = ⎢            ⎥        with  Ṙ = −[ω]ₓ R
        ⎣ −[ω]ₓR   R ⎦

Applied to `[r; v]` this gives `v_target = R·v − ω × (R·r)`, which is the
velocity relation for a rotating frame.

This is the counterpart to [`_rotation_no_rate`](@ref), and the distinction
between them is the whole point: for a spin edge the rate term is the
*dominant* velocity effect, not a correction. Using `_rotation_no_rate` here
would leave position exactly right and velocity wrong by up to `ω·r` — at
Earth's rate that is ~0.5 km/s at the surface.
"""
function _rotation_with_spin(R::AbstractMatrix, ω::Real)
    # [ω]ₓ for ω = (0, 0, ω): the cross-product matrix of a pure z rotation.
    Ω = @SMatrix [0.0  -ω   0.0
                  ω     0.0 0.0
                  0.0   0.0 0.0]
    Ṙ = -Ω * R
    return @SMatrix [R[1,1] R[1,2] R[1,3] 0.0    0.0    0.0
                     R[2,1] R[2,2] R[2,3] 0.0    0.0    0.0
                     R[3,1] R[3,2] R[3,3] 0.0    0.0    0.0
                     Ṙ[1,1] Ṙ[1,2] Ṙ[1,3] R[1,1] R[1,2] R[1,3]
                     Ṙ[2,1] Ṙ[2,2] Ṙ[2,3] R[2,1] R[2,2] R[2,3]
                     Ṙ[3,1] Ṙ[3,2] Ṙ[3,3] R[3,1] R[3,2] R[3,3]]
end

"""
    _rotation_with_rate(R::AbstractMatrix, Ṙ::AbstractMatrix) -> SMatrix{6,6,Float64,36}

Assemble the 6×6 state rotation from a rotation and its time derivative:

        ⎡ R   0 ⎤
    M = ⎢       ⎥
        ⎣ Ṙ   R ⎦

The general form. [`_rotation_no_rate`](@ref) and [`_rotation_with_spin`](@ref)
are the two special cases worth naming — `Ṙ = 0`, and `Ṙ = −[ω]ₓR` for a frame
spinning about its own z-axis. Use this where `Ṙ` is neither: an orbit-relative
frame, whose axes turn with the trajectory rather than at a constant rate.
"""
function _rotation_with_rate(R::AbstractMatrix, Ṙ::AbstractMatrix)
    return @SMatrix [R[1,1] R[1,2] R[1,3] 0.0    0.0    0.0
                     R[2,1] R[2,2] R[2,3] 0.0    0.0    0.0
                     R[3,1] R[3,2] R[3,3] 0.0    0.0    0.0
                     Ṙ[1,1] Ṙ[1,2] Ṙ[1,3] R[1,1] R[1,2] R[1,3]
                     Ṙ[2,1] Ṙ[2,2] Ṙ[2,3] R[2,1] R[2,2] R[2,3]
                     Ṙ[3,1] Ṙ[3,2] Ṙ[3,3] R[3,1] R[3,2] R[3,3]]
end
