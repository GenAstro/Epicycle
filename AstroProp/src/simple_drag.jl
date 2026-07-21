# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    Cd <: AbstractParamTag

Tag identifying the drag coefficient on a `SimpleDrag` force model.
"""
struct Cd <: AbstractParamTag end

"""
    SimpleDrag <: OrbitODE

Constant-density atmospheric drag force model.  Intended for testing and
simple mission analysis where an exponential or tabulated atmosphere is not
needed.  Keeps all tests SPICE-free and deterministic.

# Fields
- `cd::T`       — drag coefficient (dimensionless)
- `area::T`     — cross-sectional area (m²)
- `rho_const::T`— constant atmospheric density (kg/m³)

# Kernel
`SimpleDrag` uses a private kernel `_drag_accel` that is generic in `cd::T`.
`ForwardDiff` can seed `cd` with a Dual number to obtain exact ∂f/∂Cd.
The same kernel is used for ∂f/∂m (seeding `mass`).

# Notes
- `Cd` and `area` sit on `SimpleDrag` temporarily.  They will move to
  `Spacecraft` in a later sprint; only the `ModelVariable` construction line
  will change — all callers of `get_field`/`set_field!`/`accel_eval!` are
  unaffected.
"""
mutable struct SimpleDrag{T<:Real} <: OrbitODE
    cd::T
    area::T
    rho_const::T

    function SimpleDrag(; cd::Real, area::Real, rho_const::Real)
        T = promote_type(typeof(cd), typeof(area), typeof(rho_const))
        return new{T}(T(cd), T(area), T(rho_const))
    end
end

# ---------------------------------------------------------------------------
# Kernel — generic in the differentiated parameter for ForwardDiff
# ---------------------------------------------------------------------------

"""
    _drag_accel(cd, area, rho, v, mass) -> Vector

Drag acceleration kernel.  Generic in `cd` and `mass` so that
`ForwardDiff.derivative` can seed either with a Dual number for exact partials.

a_drag = -½ · cd · (area/mass) · ρ · ‖v‖ · v
"""
@inline function _drag_accel(cd::T, area, rho, v::AbstractVector, mass) where T
    return -0.5 * cd * (area / mass) * rho * norm(v) .* v
end

# ---------------------------------------------------------------------------
# accel_eval! — public force interface
# ---------------------------------------------------------------------------

function accel_eval!(m::SimpleDrag, t, y::AbstractVector, acc::AbstractVector,
                     sc::Spacecraft, params)
    v = y[4:6]
    acc[4:6] .+= _drag_accel(m.cd, m.area, m.rho_const, v, sc.mass)
    return acc
end

# ---------------------------------------------------------------------------
# state_jac! — ∂f/∂y  (analytic registration via dispatch)
# ---------------------------------------------------------------------------

function state_jac!(out::AbstractMatrix, m::SimpleDrag, t, y::AbstractVector,
                     sc::Spacecraft)
    v = y[4:6]
    # velocity block: ∂a/∂v (3×3); position block is zero
    out[4:6, 4:6] .+= ForwardDiff.jacobian(
        v_ -> _drag_accel(m.cd, m.area, m.rho_const, v_, sc.mass), v)
    return nothing
end

# ---------------------------------------------------------------------------
# param_jac! — ∂f/∂p  (analytic registration via dispatch)
# ---------------------------------------------------------------------------

function param_jac!(out::AbstractVector, m::SimpleDrag, ::Cd, t,
                     y::AbstractVector, sc::Spacecraft)
    v = y[4:6]
    out[4:6] .+= ForwardDiff.derivative(
        cd -> _drag_accel(cd, m.area, m.rho_const, v, sc.mass), m.cd)
    return nothing
end

function param_jac!(out::AbstractVector, m::SimpleDrag, ::Mass, t,
                     y::AbstractVector, sc::Spacecraft)
    v = y[4:6]
    out[4:6] .+= ForwardDiff.derivative(
        mass -> _drag_accel(m.cd, m.area, m.rho_const, v, mass), sc.mass)
    return nothing
end

# ---------------------------------------------------------------------------
# get_field / set_field! — tag interface
# ---------------------------------------------------------------------------

get_field(m::SimpleDrag, ::Cd) = m.cd

function set_field!(m::SimpleDrag, ::Cd, v::Real)
    m.cd = convert(typeof(m.cd), v)
    return nothing
end
