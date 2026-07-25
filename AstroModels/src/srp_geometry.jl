# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    AbstractSRPGeometry

Supertype for spacecraft solar-radiation-pressure geometry — the reflective shape /
coefficient data an SRP force reads from the spacecraft. Fidelity is a concrete subtype:
`CannonballSRP` now; higher-fidelity forms (e.g. an n-plate model) slot in later without
changing `Spacecraft`. Parallel to `AbstractDragGeometry`.
"""
abstract type AbstractSRPGeometry end

"""
    CannonballSRP(; c_r, srp_area)

Isotropic ("cannonball") SRP geometry — a constant reflectivity coefficient and reference
area. Stored on `Spacecraft.srp` and read by the SRP force.

# Fields
- `c_r::T`       : reflectivity coefficient (dimensionless)
- `srp_area::T`  : reference cross-sectional area [m²]
"""
struct CannonballSRP{T<:Real} <: AbstractSRPGeometry
    c_r::T
    srp_area::T
end

CannonballSRP(; c_r::Real, srp_area::Real) =
    CannonballSRP(promote(float(c_r), float(srp_area))...)

Base.show(io::IO, g::CannonballSRP) =
    print(io, "CannonballSRP(c_r = ", g.c_r, ", srp_area = ", g.srp_area, " m²)")
