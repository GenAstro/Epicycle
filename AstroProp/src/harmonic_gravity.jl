# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

# =============================================================================
# Spherical-harmonic gravity facade.
#
# Vendor-neutral user-facing wrapper around an external harmonic gravity
# implementation. Exposes only Epicycle-owned types:
#   * `HarmonicGravity(body; degree, order, coefficients=EGM96())`
#   * coefficient tags `EGM96`, …
# and forwards `accel_eval!` through an internal `ExternalForce` adapter.
# =============================================================================

using SatelliteToolboxGravityModels:
    GravityModels, IcgemFile, fetch_icgem_file
using AstroForceModels: GravityHarmonicsAstroModel

# --- Coefficient source tags ----------------------------------------------

"""Abstract supertype for a named gravity-coefficient source (e.g. EGM-96)."""
abstract type AbstractGravityCoefficients end

"""
    EGM96([path])

EGM-96 gravity coefficient source. With no argument, the ICGEM file is
fetched on demand by `SatelliteToolboxGravityModels` and cached on disk.
Supply an explicit `path` to pin a committed copy for regression tests.
"""
struct EGM96{P<:Union{Nothing,AbstractString}} <: AbstractGravityCoefficients
    path::P
end
EGM96() = EGM96(nothing)

_load_coefficients(::EGM96{Nothing}) =
    GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
_load_coefficients(c::EGM96{<:AbstractString}) =
    GravityModels.load(IcgemFile, c.path)

# --- Public facade --------------------------------------------------------

"""
    HarmonicGravity(body; degree, order, coefficients=EGM96())

Spherical-harmonic gravity force model for `body` (currently Earth-tested).
`degree=0, order=0` reduces to a pure central-body monopole and is useful as
a plumbing sanity check against `PointMassGravity`.

The active `AstroUniverse.eop()` table is captured at construction time, so
EOP must be configured before building the model. Internally this composes
a HAMMERHEAD `AstroForceModels.GravityHarmonicsAstroModel`; that type is
intentionally not part of Epicycle's public API.
"""
struct HarmonicGravity{F} <: AbstractGravityForce
    central_body::CelestialBody
    degree::Int
    order::Int
    ext::ExternalForce{F}
end

function HarmonicGravity(body::CelestialBody;
                         degree::Integer,
                         order::Integer,
                         coefficients::AbstractGravityCoefficients = EGM96())
    grav_model = _load_coefficients(coefficients)
    hh = GravityHarmonicsAstroModel(;
        gravity_model = grav_model,
        eop_data      = AstroUniverse.eop(),
        degree        = Int(degree),
        order         = Int(order),
    )
    return HarmonicGravity(body, Int(degree), Int(order),
                           ExternalForce(hh, true))   # primary gravity
end

@inline accel_eval!(g::HarmonicGravity, t::Time, y::AbstractVector,
                    dy::AbstractVector, sc::Spacecraft, params) =
    accel_eval!(g.ext, t, y, dy, sc, params)
