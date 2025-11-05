# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

using LinearAlgebra

"""
    cart_to_sphazfpa(cart::Vector{<:Real}; tol::Float64 = 1e-12) 

Convert a Cartesian state to Spherical AZ-FPA representation.

# Arguments
- `cart::Vector{<:Real}`: Cartesian state `[x, y, z, vx, vy, vz]`
- `tol::Float64`: Numerical tolerance for singularity checks (default: `1e-12`)

# Returns
A 6-element Spherical AZ-FPA state `[r, λ, δ, v, αₚ, ψ]`:
- `r`   : radial distance [length]
- `λ`   : right ascension [rad]
- `δ`   : declination [rad]
- `v`   : velocity magnitude [length/time]
- `αₚ`  : flight path azimuth [rad]
- `ψ`   : flight path angle [rad]

# Notes
- Returns `fill(NaN, 6)` if `r` or `v` are near zero or orbit is singular
- All angles are in radians.

# Examples
```julia
cart = [6778.0, 0.0, 0.0, 0.0, 7.66, 0.0]
sphazfpa = cart_to_sphazfpa(cart)
```
"""
function cart_to_sphazfpa(cart::Vector{<:Real}; tol::Float64 = 1e-12)
    if length(cart) != 6
        error("Input vector must have six elements: [x, y, z, vx, vy, vz]")
    end

    r̄ = cart[1:3]
    v̄ = cart[4:6]

    r = norm(r̄)
    if r < tol
        @warn "Conversion failed: Position magnitude r = $r is below tolerance."
        return fill(NaN, 6)
    end

    v = norm(v̄)
    if v < tol
        @warn "Conversion failed: Velocity magnitude v = $v is below tolerance."
        return fill(NaN, 6)
    end

    λ = atan(r̄[2], r̄[1])
    δ = asin(r̄[3] / r)
    ψ = acos(clamp(dot(r̄, v̄) / (r * v), -1.0, 1.0))  

    # Build local spherical frame
    x̂ = [cos(δ) * cos(λ),  cos(δ) * sin(λ), sin(δ)]                      # radial
    ŷ = [cos(λ + π/2),     sin(λ + π/2),     0.0]                        # east
    ẑ = [-sin(δ) * cos(λ), -sin(δ) * sin(λ), cos(δ)]                     # north

    R_li = hcat(x̂, ŷ, ẑ)'  # transpose → from inertial to local frame

    v_local = R_li * v̄
    αₚ = atan(v_local[2], v_local[3])  # azimuth from north toward east

    return [r, λ, δ, v, αₚ, ψ]
end
