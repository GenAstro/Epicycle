# Copyright (C) 2025 Gen Astro LLC
SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    cart_to_mee(cart::Vector{<:Real}, μ::Real; j::Float64 = 1.0)

Convert Cartesian state to Modified Equinoctial Elements (MEE).

# Arguments
- `cart::Vector{<:Real}`: 6-element vector `[x, y, z, vx, vy, vz]`
- `μ::Real`: Gravitational parameter
- `j::Float64=1.0`: Optional constant (1 for prograde, -1 for retrograde), defaults to `1.0`
- `tol::Float64`: Optional tolerance for singularity checking

# Returns
- A 6-element vector `[p, f, g, h, k, L]` representing the modified equinoctial elements.

# Examples
```julia
cart = [6778.0, 0.0, 0.0, 0.0, 7.66, 0.0]
mee = cart_to_mee(cart, 398600.4418)
```
"""
function cart_to_mee(cart::Vector{<:Real}, μ::Real; j::Float64 = 1.0, tol::Float64 = 1e-12)
    if length(cart) != 6
        error("Input vector must have exactly six elements: [x, y, z, vx, vy, vz].")
    end

    # Validate j
    if j ∉ (-1.0, 1.0)
        error("Invalid value for j: must be 1.0 or -1.0")
    end

    # Split input vector into position and velocity
    r̄ = cart[1:3]
    v̄ = cart[4:6]
    r = norm(r̄)
   
    # Angular momentum vector and magnitude
    h̄ = cross(r̄, v̄)
    h = norm(h̄)

    # Unit vectors
    r̂ = r == 0 ? zeros(3) : r̄ / r
    if h == 0
        ĥ = zeros(3)
        v_hat = zeros(3)
    else
        ĥ = h̄ / h
        v_hat = (r * v̄ - dot(r̄, v̄) * r̄ / r) / h
    end

    # Eccentricity vector
    ē = cross(v̄, h̄) / μ - r̂

    # Semi-latus rectum
    p = h^2 / μ
    if p < 0
        error("Semi-latus rectum must be greater than 0")
    end

    # Avoid singularity when computing h and k
    denom = 1.0 + ĥ[3] * j
    if abs(denom) < tol
        @warn "Singularity computing h and k while computing mee elements"
        return fill(NaN, 6)
    end

    # Reference frame for computing f, g, h, k
    fx = 1 - ĥ[1]^2 / denom
    fy = -ĥ[1] * ĥ[2] / denom
    fz = -ĥ[1] * j
    f̂ = [fx, fy, fz]
    ĝ = cross(ĥ, f̂)

    # Compute modified equinoctial elements
    f = dot(ē, f̂)
    g = dot(ē, ĝ)
    h = -ĥ[2] / denom
    k =  ĥ[1] / denom

    # Compute true longitude L
    sinl = r̂[2] - v_hat[1]
    cosl = r̂[1] + v_hat[2]
    L = mod(atan(sinl, cosl), 2π)

    return [p, f, g, h, k, L]
end
