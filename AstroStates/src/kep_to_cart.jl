# Copyright 2025 Gen Astro LLC. All Rights Reserved.
#
# This software is licensed under the GNU AGPL v3.0,
# WITHOUT ANY WARRANTY, including implied warranties of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# This file may also be used under a commercial license,
# if one has been purchased from Gen Astro LLC.
#
# By modifying this software, you agree to the terms of the
# Gen Astro LLC Contributor License Agreement.

using LinearAlgebra

"""
    kep_to_cart(state::Vector{<:Real}, μ::Real; tol::Float64=1e-12)

Convert a Keplerian state vector to a Cartesian state vector.

# Arguments
- `state::Vector{<:Real}`: Keplerian elements `[a, e, i, Ω, ω, ν]`
- `μ`: Gravitational parameter
- `tol`: Tolerance for singularities like p ≈ 0 (default: 1e-12)
- `a`: semi-major axis
- `e`: eccentricity
- `i`: inclination
- `Ω`: right ascension of ascending node
- `ω`: argument of periapsis
- `ν`: true anomaly

# Returns
A 6-element vector `[x, y, z, vx, vy, vz]` representing Cartesian position and velocity.

# Examples
```julia
kep = [7000.0, 0.01, π/4, 0.0, 0.0, π/3]
cart = kep_to_cart(kep, 398600.4418)
```

# Notes
- Angles must be in radians.
- Dimensional quantities must be consistent units with μ.
- Returns a vector of `NaN`s if conversion is undefined.
"""
function kep_to_cart(state::Vector{<:Real}, μ::Real; tol::Float64=1e-12)
    if length(state) != 6
        error("Input vector must have exactly six elements: a, e, i, Ω, ω, ν.")
    end

    if μ < tol
        @warn "Conversion Failed: μ < tolerance."
        return fill(NaN, 6)
    end

    # Unpack the elements
    a, e, i, Ω, ω, ν = state

    # Compute semi-latus rectum: p = a * (1 - e²)
    p = a * (1.0 - e^2)

    # Check for degenerate orbit (e.g., parabolic or collapsed)
    if p < tol || abs(1-e) < tol
        @warn "Conversion Failed: Orbit is parabolic or singular."
        return fill(NaN, 6)
    end

    # Compute radial distance: r = p / (1 + e * cos(ν))
    r = p / (1.0 + e * cos(ν))

    # Position and velocity in perifocal frame 
    factor = sqrt(μ / p)
    r̄ₚ = [r * cos(ν), r * sin(ν), 0.0]
    v̄ₚ = [-factor * sin(ν), factor * (e + cos(ν)), 0.0]

    # Precompute sines and cosines for rotation matrix
    cos_Ω, sin_Ω = cos(Ω), sin(Ω)
    cos_ω, sin_ω = cos(ω), sin(ω)
    cos_i, sin_i = cos(i), sin(i)

    # Rotation matrix from perifocal to inertial
    R = [
        cos_ω * cos_Ω - sin_ω * cos_i * sin_Ω   -sin_ω * cos_Ω - cos_ω * cos_i * sin_Ω   sin_i * sin_Ω;
        cos_ω * sin_Ω + sin_ω * cos_i * cos_Ω   -sin_ω * sin_Ω + cos_ω * cos_i * cos_Ω  -sin_i * cos_Ω;
        sin_ω * sin_i                                    cos_ω * sin_i                   cos_i
    ]

    # Rotate position and velocity from perifocal to inertial frame
    pos = R * r̄ₚ
    vel = R * v̄ₚ 

    return vcat(pos, vel)
end
