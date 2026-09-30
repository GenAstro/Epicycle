# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using LinearAlgebra

"""
    sphazfpa_to_cart(spherical::AbstractVector{<:Real}) -> Vector

Convert a Spherical AZ-FPA state to a Cartesian state vector.

# Arguments
- `spherical`: Spherical AZ-FPA state vector `[r, λ, δ, v, αₚ, ψ]`
    - `r`   : radial distance [length]
    - `λ`   : right ascension [rad]
    - `δ`   : declination [rad]
    - `v`   : velocity magnitude [length/time]
    - `αₚ`  : flight path azimuth (angle east of north in the local horizontal plane) [rad]
    - `ψ`   : flight path angle, measured from the radial direction [rad]; π/2 is horizontal
              flight. This is GMAT's convention.

# Returns
A 6-element Cartesian state vector `[x, y, z, vx, vy, vz]`.

# Notes
- All angles must be in radians.

# Examples
```julia
sphazfpa = [6478.0, 0.0, π/4, 7.5, π/4, π/2]   # horizontal flight toward the northeast
cart = sphazfpa_to_cart(sphazfpa)
```
"""
function sphazfpa_to_cart(spherical::AbstractVector{<:Real})
    if length(spherical) != 6
        error("Input vector must have six elements: [r, λ, δ, v, αₚ, ψ]")
    end
    T = float(eltype(spherical))

    r, λ, δ, v, αₚ, ψ = spherical

    # Precompute trigonometric terms
    sinδ, cosδ = sincos(δ)
    sinλ, cosλ = sincos(λ)
    sinψ, cosψ = sincos(ψ)
    sinα, cosα = sincos(αₚ)

    # Position components
    x = r * cosδ * cosλ
    y = r * cosδ * sinλ
    z = r * sinδ

    # Velocity: cos ψ along the radial direction, sin ψ in the horizontal plane at azimuth αₚ
    vx = v * ( cosψ * cosδ * cosλ -
               sinψ * (sinα * sinλ + cosα * sinδ * cosλ) )
    vy = v * ( cosψ * cosδ * sinλ +
               sinψ * (sinα * cosλ - cosα * sinδ * sinλ) )
    vz = v * ( cosψ * sinδ + sinψ * cosα * cosδ )

    return T[x, y, z, vx, vy, vz]
end
