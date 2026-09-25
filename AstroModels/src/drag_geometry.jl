# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

"""
    AbstractDragGeometry

`AbstractDragGeometry` is the common type for the shape and coefficient data a
drag force reads from a spacecraft. `SphericalDrag` provides the isotropic model;
other fidelity levels use additional concrete subtypes.

# Example
```jldoctest
SphericalDrag(c_d = 2.2, drag_area = 4.0) isa AbstractDragGeometry

# output
true
```
"""
abstract type AbstractDragGeometry end

"""
    SphericalDrag(; c_d, drag_area)

`SphericalDrag` defines isotropic drag geometry with a constant drag coefficient
and reference area. A spacecraft stores the model in `drag`.

# Fields
- `c_d::T`       : drag coefficient (dimensionless)
- `drag_area::T` : reference cross-sectional area [m²]

# Example
```julia
drag = SphericalDrag(c_d = 2.2, drag_area = 4.0)
```
"""
struct SphericalDrag{T<:Real} <: AbstractDragGeometry
    c_d::T
    drag_area::T
end

SphericalDrag(; c_d::Real, drag_area::Real) =
    SphericalDrag(promote(float(c_d), float(drag_area))...)

Base.show(io::IO, g::SphericalDrag) =
    print(io, "SphericalDrag(c_d = ", g.c_d, ", drag_area = ", g.drag_area, " m²)")

function Base.show(io::IO, ::MIME"text/plain", g::SphericalDrag)
    println(io, "SphericalDrag:")
    println(io, "  c_d  = ", g.c_d)
    print(io,   "  area = ", g.drag_area, " m²")
end
