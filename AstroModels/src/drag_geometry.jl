# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    AbstractDragGeometry

Supertype for spacecraft drag geometry — the physical shape/coefficient data that a
drag force reads from the spacecraft. Fidelity is a concrete subtype: `CannonballDrag`
now; higher-fidelity forms (e.g. an n-plate model loaded from a file) slot in later
without changing `Spacecraft`.
"""
abstract type AbstractDragGeometry end

"""
    CannonballDrag(; c_d, drag_area)

Isotropic ("cannonball") drag geometry — a constant drag coefficient and reference
area. Stored on `Spacecraft.drag` and read by the drag force.

# Fields
- `c_d::T`       : drag coefficient (dimensionless)
- `drag_area::T` : reference cross-sectional area [m²]
"""
struct CannonballDrag{T<:Real} <: AbstractDragGeometry
    c_d::T
    drag_area::T
end

CannonballDrag(; c_d::Real, drag_area::Real) =
    CannonballDrag(promote(float(c_d), float(drag_area))...)

Base.show(io::IO, g::CannonballDrag) =
    print(io, "CannonballDrag(c_d = ", g.c_d, ", drag_area = ", g.drag_area, " m²)")
