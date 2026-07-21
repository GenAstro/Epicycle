# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    Cr <: AbstractParamTag

Tag identifying the radiation pressure coefficient on a `SimpleSRP` force model.
"""
struct Cr <: AbstractParamTag end

"""
    SimpleSRP <: OrbitODE

Simple solar radiation pressure force model with a fixed sun direction vector.
Intended for testing and simple mission analysis — no ephemeris required, so
all tests remain deterministic and SPICE-free.

# Fields
- `cr::T`               — radiation pressure coefficient (dimensionless)
- `area::T`             — cross-sectional area (m²)
- `P_srp::T`            — solar radiation pressure at 1 AU (N/m²)
- `sun_dir::Vector{T}`  — fixed unit vector toward the sun (inertial frame)

# Kernel
`SimpleSRP` uses a private kernel `_srp_accel` that is generic in `cr::T`.
`ForwardDiff` can seed `cr` or `mass` with a Dual number for exact partials.
"""
mutable struct SimpleSRP{T<:Real} <: OrbitODE
    cr::T
    area::T
    P_srp::T
    sun_dir::Vector{T}

    function SimpleSRP(; cr::Real, area::Real, P_srp::Real,
                         sun_dir::AbstractVector)
        T = promote_type(typeof(cr), typeof(area), typeof(P_srp), eltype(sun_dir))
        return new{T}(T(cr), T(area), T(P_srp), Vector{T}(sun_dir))
    end
end

# ---------------------------------------------------------------------------
# Kernel — generic in the differentiated parameter for ForwardDiff
# ---------------------------------------------------------------------------

"""
    _srp_accel(cr, area, P_srp, sun_dir, mass) -> Vector

SRP acceleration kernel.  Generic in `cr` and `mass` so that
`ForwardDiff.derivative` can seed either for exact partials.

a_srp = (cr · area · P_srp / mass) · sun_dir
"""
@inline function _srp_accel(cr::T, area, P_srp, sun_dir::AbstractVector,
                              mass) where T
    return (cr * area * P_srp / mass) .* sun_dir
end

# ---------------------------------------------------------------------------
# accel_eval! — public force interface
# ---------------------------------------------------------------------------

function accel_eval!(m::SimpleSRP, t, y::AbstractVector, acc::AbstractVector,
                     sc::Spacecraft, params)
    acc[4:6] .+= _srp_accel(m.cr, m.area, m.P_srp, m.sun_dir, sc.mass)
    return acc
end

# ---------------------------------------------------------------------------
# state_jac! — ∂f/∂y  (analytic registration via dispatch)
# SRP is independent of state — register an explicit no-op to avoid the
# ForwardDiff.jacobian fallback computing zeros expensively.
# ---------------------------------------------------------------------------

function state_jac!(out::AbstractMatrix, m::SimpleSRP, t, y::AbstractVector,
                     sc::Spacecraft)
    return nothing  # SRP is state-independent
end

# ---------------------------------------------------------------------------
# param_jac! — ∂f/∂p  (analytic registration via dispatch)
# ---------------------------------------------------------------------------

function param_jac!(out::AbstractVector, m::SimpleSRP, ::Cr, t,
                     y::AbstractVector, sc::Spacecraft)
    out[4:6] .+= ForwardDiff.derivative(
        cr -> _srp_accel(cr, m.area, m.P_srp, m.sun_dir, sc.mass), m.cr)
    return nothing
end

function param_jac!(out::AbstractVector, m::SimpleSRP, ::Mass, t,
                     y::AbstractVector, sc::Spacecraft)
    out[4:6] .+= ForwardDiff.derivative(
        mass -> _srp_accel(m.cr, m.area, m.P_srp, m.sun_dir, mass), sc.mass)
    return nothing
end

# ---------------------------------------------------------------------------
# get_field / set_field! — tag interface
# ---------------------------------------------------------------------------

get_field(m::SimpleSRP, ::Cr) = m.cr

function set_field!(m::SimpleSRP, ::Cr, v::Real)
    m.cr = convert(typeof(m.cr), v)
    return nothing
end
