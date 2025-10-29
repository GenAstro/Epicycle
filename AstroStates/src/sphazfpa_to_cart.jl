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
    sphazfpa_to_cart(spherical::Vector{<:Real}) -> Vector{<:Real}

Convert a Spherical AZ-FPA state to a Cartesian state vector.

# Arguments
- `spherical::Vector{<:Real}`: Spherical AZ-FPA state vector `[r, λ, δ, v, αₚ, ψ]`
    - `r`   : radial distance [length]
    - `λ`   : right ascension [rad]
    - `δ`   : declination [rad]
    - `v`   : velocity magnitude [length/time]
    - `αₚ`  : flight path azimuth (angle east of north in local horizon) [rad]
    - `ψ`   : flight path angle (angle above local horizon) [rad]

# Returns
A 6-element Cartesian state vector `[x, y, z, vx, vy, vz]`.

# Notes
- All angles must be in radians.
- Velocity frame uses local vertical/local horizontal.

# Examples
```julia
sphazfpa = [6478.0, 0.0, π/4, 7.5, π/4, π/4]
cart = sphazfpa_to_cart(sphazfpa)
```
"""
function sphazfpa_to_cart(spherical::Vector{<:Real})
    if length(spherical) != 6
        error("Input vector must have six elements: [r, λ, δ, v, αₚ, ψ]")
    end

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

    # Velocity in components
    vx = v * ( cosψ * cosδ * cosλ -
               sinψ * (sinα * sinλ + cosα * sinδ * cosλ) )
    vy = v * ( cosψ * cosδ * sinλ +
               sinψ * (sinα * cosλ - cosα * sinδ * sinλ) )
    vz = v * ( cosψ * sinδ + sinψ * cosα * cosδ )

    return [x, y, z, vx, vy, vz]
end
