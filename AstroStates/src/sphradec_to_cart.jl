# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

using LinearAlgebra

"""
    sphradec_to_cart(state::Vector{<:Real}) -> Vector{<:Real}

Convert a spherical RA/Dec state `[r, λᵣ, δᵣ, v, λᵥ, δᵥ]` to a Cartesian state
`[x, y, z, vx, vy, vz]`.

Arguments
- state::Vector{<:Real}: length-6 vector where:
  - r   = position magnitude
  - λᵣ  = right ascension of position (radians)
  - δᵣ  = declination of position (radians)
  - v   = velocity magnitude
  - λᵥ  = azimuth of velocity direction (radians)
  - δᵥ  = elevation of velocity direction (radians)

Returns
- Vector{<:Real}: length-6 Cartesian state `[x, y, z, vx, vy, vz]`.

Notes
- All angles are in radians.
- Units must be consistent between position and velocity components.
- Right ascension typically in [0, 2π); declination in [-π/2, π/2].

Examples
```julia
sphradec = [7000.0, 0.0, 0.0, 7.5, π/2, 0.0]
cart = sphradec_to_cart(sphradec)
```
"""
function sphradec_to_cart(state::Vector{<:Real})
    if length(state) != 6
        error("Input vector must have exactly six elements: [r, λᵣ, δᵣ, v, λᵥ, δᵥ].")
    end

    r, λᵣ, δᵣ, v, λᵥ, δᵥ = state

    # Position conversion: spherical to Cartesian
    x = r * cos(δᵣ) * cos(λᵣ)
    y = r * cos(δᵣ) * sin(λᵣ)
    z = r * sin(δᵣ)

    # Velocity conversion: spherical to Cartesian
    vx = v * cos(λᵥ) * cos(δᵥ)
    vy = v * sin(λᵥ) * cos(δᵥ)
    vz = v * sin(δᵥ)  

    return [x, y, z, vx, vy, vz]
end
