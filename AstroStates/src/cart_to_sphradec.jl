# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using LinearAlgebra

"""
    cart_to_sphradec(state::AbstractVector{<:Real}; tol::Real=1e-12)

Convert a Cartesian state vector to a spherical RA/DEC state vector.

# Arguments
- `state`: A 6-element vector `[x, y, z, vx, vy, vz]` representing Cartesian position and velocity.

# Returns
A 6-element vector `[r, ra, dec, v, rav, decv]`, the field order of `SphericalRADECState`:
- `r`  : magnitude of position vector
- `λᵣ` : right ascension (radians), in [0, 2π)
- `δᵣ` : declination (radians), in [-π/2, π/2]
- `v`  : magnitude of velocity
- `λᵥ` : right ascension of the velocity (radians), in [0, 2π)
- `δᵥ` : declination of the velocity (radians), in [-π/2, π/2]

A zero position or velocity logs a warning and returns `NaN`s.

# Notes
- Assumes all angles are in radians.
- Units must be consistent between position and velocity components.

# Examples
```julia
cart = [6778.0, 0.0, 0.0, 0.0, 7.66, 0.0]
sphradec = cart_to_sphradec(cart)
```
"""
function cart_to_sphradec(state::AbstractVector{<:Real}; tol::Real=1e-12)
    if length(state) != 6
        error("Input vector must have six elements: [x, y, z, vx, vy, vz].")
    end
    T = float(eltype(state))

    # Unpack Cartesian position and velocity
    x, y, z, vx, vy, vz = state

    # Magnitude of position
    r = sqrt(x^2 + y^2 + z^2)
    if r < tol
        @warn "Conversion failed: Radius is zero."
        return fill(T(NaN), 6)
    end

    # Spherical angles of the position
    λᵣ = _wrap_2pi(atan(y, x))
    δᵣ = atan(z, hypot(x, y))

    # Magnitude of velocity
    v = sqrt(vx^2 + vy^2 + vz^2)
    if v < tol
        @warn "Conversion failed: Velocity is zero."
        return fill(T(NaN), 6)
    end

    # Spherical angles of the velocity
    λᵥ = _wrap_2pi(atan(vy, vx))
    δᵥ = atan(vz, hypot(vx, vy))

    return T[r, λᵣ, δᵣ, v, λᵥ, δᵥ]
end
