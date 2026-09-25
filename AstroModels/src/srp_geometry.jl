# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

"""
    AbstractSRPGeometry

`AbstractSRPGeometry` is the common type for the reflective shape and coefficient
data a solar-radiation-pressure force reads from a spacecraft. `SphericalSRP`
provides the isotropic model; other fidelity levels use additional concrete subtypes.

# Example
```jldoctest
SphericalSRP(c_r = 1.8, srp_area = 10.0) isa AbstractSRPGeometry

# output
true
```
"""
abstract type AbstractSRPGeometry end

"""
    SphericalSRP(; c_r, srp_area)

`SphericalSRP` defines isotropic solar-radiation-pressure geometry with a
constant reflectivity coefficient and reference area. A spacecraft stores the
model in `srp`.

# Fields
- `c_r::T`       : reflectivity coefficient (dimensionless)
- `srp_area::T`  : reference cross-sectional area [m²]

# Example
```julia
srp = SphericalSRP(c_r = 1.8, srp_area = 10.0)
```
"""
struct SphericalSRP{T<:Real} <: AbstractSRPGeometry
    c_r::T
    srp_area::T
end

SphericalSRP(; c_r::Real, srp_area::Real) =
    SphericalSRP(promote(float(c_r), float(srp_area))...)

Base.show(io::IO, g::SphericalSRP) =
    print(io, "SphericalSRP(c_r = ", g.c_r, ", srp_area = ", g.srp_area, " m²)")

function Base.show(io::IO, ::MIME"text/plain", g::SphericalSRP)
    println(io, "SphericalSRP:")
    println(io, "  c_r  = ", g.c_r)
    print(io,   "  area = ", g.srp_area, " m²")
end
