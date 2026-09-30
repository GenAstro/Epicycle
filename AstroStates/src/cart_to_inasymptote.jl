# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    cart_to_inasymptote(cart::AbstractVector{<:Real}, μ::Real; tol::Real = 1e-12) -> Vector

Convert a Cartesian state to incoming asymptote parameters.

# Arguments
- `cart`: Cartesian state `[x, y, z, vx, vy, vz]`
- `μ::Real`: gravitational parameter
- `tol::Real`: tolerance for detecting singularities

# Returns
A 6-element vector: `[rₚ, C₃, λₐ, δₐ, θᵦ, ν]`
- rₚ : radius of periapsis
- C₃ : characteristic energy
- λₐ : right ascension of incoming asymptote, in [0, 2π)
- δₐ : declination of incoming asymptote, in [-π/2, π/2]
- θᵦ : B-plane angle, in [0, 2π)
- ν : true anomaly, in [0, 2π)

For an elliptic orbit (C₃ < 0) there is no asymptote, and the apoapsis direction is used in its
place. A circular or parabolic orbit, or an asymptote along the z-axis, logs a warning and
returns `NaN`s.

# Notes
- Angles are in radians.
- Dimensional quantities are consistent units with μ.

# Examples
```julia
cart = [10000.0, 0.0, 0.0, 0.0, 12.0, 0.0]  # Hyperbolic trajectory
inasym = cart_to_inasymptote(cart, 398600.4418)
```
"""
cart_to_inasymptote(cart::AbstractVector{<:Real}, μ::Real; tol::Real = 1e-12) =
    _cart_to_asymptote(cart, μ, -1, tol)
