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
    cart_to_sphradec(state::Vector{<:Real}; tol::Float64=1e-12) 

Convert a Cartesian state vector to a spherical RA/DEC state vector.

# Arguments
- `state::Vector{<:Real}`: A 6-element vector `[x, y, z, vx, vy, vz]` representing Cartesian position and velocity.

# Returns
A 6-element vector `[r, ra, dec, v, vra, vdec]` where:
- `r`  : magnitude of position vector 
- `λᵣ` : right ascension (radians)
- `δᵣ` : declination (radians)
- `v`  : magnitude of velocity
- `λᵥ` : azimuthal direction of velocity (radians)
- `δᵥ` : elevation direction of velocity (radians)

# Notes
- Assumes all angles are in radians.
- Units must be consistent between position and velocity components.

# Examples
```julia
cart = [6778.0, 0.0, 0.0, 0.0, 7.66, 0.0]
sphradec = cart_to_sphradec(cart)
```
"""
function cart_to_sphradec(state::Vector{<:Real}; tol::Float64=1e-12)
    if length(state) != 6
        error("Input vector must have six elements: [x, y, z, vx, vy, vz].")
    end

    # Unpack Cartesian position and velocity
    x, y, z, vx, vy, vz = state

    # Compute magnitude of position
    r = sqrt(x^2 + y^2 + z^2)
    if r < tol
        @warn "Conversion failed: Radius is zero."
        return fill(NaN, 6)
    end

    # Compute spherical angles from position
    λᵣ = atan(y, x)  
    δᵣ = asin(z / r) 

    # Compute magnitude of velocity
    v = sqrt(vx^2 + vy^2 + vz^2)
    if v < tol
        @warn "Conversion failed: Velocity is zero."
        return fill(NaN, 6)
    end

    # Compute spherical angles from velocity
    λᵥ = atan(vy, vx)  
    δᵥ = asin(vz / v)

    return [r, λᵣ, δᵣ, v, λᵥ, δᵥ]
end
