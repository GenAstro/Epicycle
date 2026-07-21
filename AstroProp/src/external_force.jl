# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

# =============================================================================
# External force-catalog adapter (internal).
#
# Thin glue between AstroProp's force interface (`accel_eval!`) and external
# force catalogs that expose an `acceleration(u, p, t, model)` API in J2000
# inertial km / km·s⁻². Currently consumes HAMMERHEAD `AstroForceModels.jl`
# models, but nothing in this file mentions that vendor by type — the public
# Epicycle API only ever sees `HarmonicGravity`, `CannonballDrag`, etc.
# =============================================================================

using AstroForceModels: acceleration, AbstractAstroForceModel
using ComponentArrays: ComponentVector
using StaticArrays: SVector

"""
    AbstractGravityForce <: OrbitODE

Common supertype for any gravity force model recognised by `ForceModel`'s
central-body inference (`_find_center`). Subtypes must expose a
`central_body::CelestialBody` field.
"""
abstract type AbstractGravityForce <: OrbitODE end

"""
    ExternalForce{F} <: OrbitODE

Internal adapter that lifts a vendor force model (subtype of
`AstroForceModels.AbstractAstroForceModel`) into AstroProp's `accel_eval!`
contract. Not exported. Wrapped by user-facing facades such as
`HarmonicGravity`.

`writes_kinematic` controls whether `accel_eval!` writes `dy[1:3] = v`. Set
to `true` for the primary gravity model and `false` for additive
perturbations (drag, SRP, third-body) so they don't clobber the kinematic
row.
"""
struct ExternalForce{F<:AbstractAstroForceModel} <: OrbitODE
    model::F
    writes_kinematic::Bool
end

ExternalForce(model::AbstractAstroForceModel) = ExternalForce(model, false)

"""
    _hh_params(t::Time) -> ComponentVector

Build the minimal parameter bag the HAMMERHEAD acceleration kernels expect.
They access `p.JD` and form `current_jd(p, t_rel) = p.JD + t_rel/86400`; by
pinning `t_rel = 0` and `p.JD = t_utc`, we sidestep any drift between
AstroProp's `Time` and the catalog's epoch.
"""
@inline _hh_params(t::Time) = ComponentVector(JD = t.utc.jd)

function accel_eval!(m::ExternalForce, t::Time, y::AbstractVector,
                     dy::AbstractVector, sc::Spacecraft, params)
    p_hh = _hh_params(t)
    a    = acceleration(SVector{6}(y[1], y[2], y[3], y[4], y[5], y[6]),
                        p_hh, 0.0, m.model)
    if m.writes_kinematic
        dy[1] = y[4]; dy[2] = y[5]; dy[3] = y[6]
    end
    dy[4] += a[1]; dy[5] += a[2]; dy[6] += a[3]
    return dy
end
