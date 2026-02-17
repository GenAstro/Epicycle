# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

using LinearAlgebra

"""
    cart_to_equinoctial(cart::Vector{<:Real}, μ::Real; tol::Float64 = 1e-12)

Convert a Cartesian state vector to equinoctial orbital elements.

# Arguments
- `cart::Vector{<:Real}`: Cartesian state `[x, y, z, vx, vy, vz]`
- `μ::Real`: gravitational parameter [length³/time²]
- `tol::Float64`: tolerance for singularity detection (default = 1e-12)

# Returns
Equinoctial state `[a, h, k, p, q, λ]`:
- `a` : semi-major axis [length]
- `h` : e⋅g component of eccentricity vector
- `k` : e⋅f component of eccentricity vector
- `p` : tan(i/2)⋅cos(Ω)
- `q` : tan(i/2)⋅sin(Ω)
- `λ` : mean longitude [rad]

# Notes
- Fails gracefully with `NaN` if orbit is hyperbolic, singular, or edge-case unsupported.
- All angles are in radians. Units consistent with `μ`.
- Note that in most cases in states, h is the magnitude of angular momentum.  But, not
  for equinoctial elements.  

# Examples
```julia
cart = [6778.0, 0.0, 0.0, 0.0, 7.66, 0.0]
equinoctial = cart_to_equinoctial(cart, 398600.4418)
```
"""
function cart_to_equinoctial(cart::Vector{<:Real}, μ::Real; tol::Float64 = 1e-12)
    if length(cart) != 6
        error("Input must be a 6-element Cartesian state vector [x, y, z, vx, vy, vz] .")
    end

    r̄ = cart[1:3]
    v̄ = cart[4:6]

    r = norm(r̄)
    v = norm(v̄)

    if r < tol
        @warn "Conversion failed: Position magnitude r = $r less than tol."
        return fill(NaN, 6)
    end
    if μ < tol
        @warn "Conversion failed: Gravitational parameter μ = $μ less than tol."
        return fill(NaN, 6)
    end

    # Angular momentum vector and check for radial/degenerate orbit
    ang_mom_vec = cross(r̄, v̄)
    h_mag = norm(ang_mom_vec)
    
    if h_mag < tol
        @warn "Conversion failed: Angular momentum near zero (radial or degenerate orbit)."
        return fill(NaN, 6)
    end

    # Eccentricity vector and magnitude
    ē = ((v^2 - μ / r) * r̄ - dot(r̄, v̄) * v̄) / μ
    e = norm(ē)

    # Specific energy and semi-major axis
    ξ = v^2 / 2 - μ / r
    a = -μ / (2 * ξ)
    
    if a < tol
        @warn "Conversion failed: Orbit is parabolic or hyperbolic (a = $a)."
        return fill(NaN, 6)
    end

    # Angular momentum unit vector and inclination
    unit_ang_mom = ang_mom_vec / h_mag
    i = acos(clamp(ang_mom_vec[3] / h_mag, -1.0, 1.0))

    if abs(i - π) < tol
        @warn "Conversion failed: Equinoctial elements not defined for i ≈ π."
        return fill(NaN, 6)
    end

    # Equinoctial reference vectors f and g
    j = 1  # For general case; would be -1 if i = π, but not supported
    denom = 1 + unit_ang_mom[3]^j
    f̂ = [
        1 - (unit_ang_mom[1]^2) / denom,
        -unit_ang_mom[1] * unit_ang_mom[2] / denom,
        -unit_ang_mom[1]^j
    ]
    f̂ = f̂/norm(f̂)
    ĝ = cross(unit_ang_mom, f̂)
    ĝ = ĝ/norm(ĝ)

    # Project eccentricity vector onto equinoctial frame
    h = dot(ē, ĝ)
    k = dot(ē, f̂)
    p = unit_ang_mom[1] / denom
    q = -unit_ang_mom[2] / denom

    # Compute true longitude
    X1 = dot(r̄, f̂)
    Y1 = dot(r̄, ĝ)
    sqrt1 = sqrt(1 - h^2 - k^2)
    β = 1 / (1 + sqrt1)

    cosF = k + ((1 - k^2 * β) * X1 - h * k * β * Y1) / (a * sqrt1)
    sinF = h + ((1 - h^2 * β) * Y1 - h * k * β * X1) / (a * sqrt1)
    F = atan(sinF, cosF)
    F = F < 0 ? F + 2π : F  

    λ = F + h * cosF - k * sinF

    return [a, h, k, p, q, λ]
end
