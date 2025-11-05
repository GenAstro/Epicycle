# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

"""
    function cart_to_kep(cart::Vector{<:Real}, μ::Real; tol::Float64=1e-12)

Convert a Cartesian state vector to Keplerian orbital elements.

# Arguments
- `cart::Vector{Float64}`: Cartesian state `[x, y, z, vx, vy, vz]`
- `μ::Real`: Gravitational parameter
- `tol::Float64`: Tolerance for detecting circular or equatorial orbits (default: `1e-12`)

# Returns
A vector `[a, e, i, Ω, ω, ν]` where:
- `a`: semi-major axis
- `e`: eccentricity
- `i`: inclination
- `Ω`: right ascension of ascending node (RAAN)
- `ω`: argument of periapsis
- `ν`: true anomaly

# Notes
- Angles must be in radians.
- Dimensional quantities must use consistent units with μ.

# Examples
```julia
cart = [6778.0, 0.0, 0.0, 0.0, 7.66, 0.0]
kep = cart_to_kep(cart, 398600.4418)
```
"""
function cart_to_kep(cart::Vector{<:Real}, μ::Real; tol::Float64=1e-12)
    if length(cart) != 6
        error("Input vector must have exactly six elements: [x, y, z, vx, vy, vz].")
    end

    if μ < tol
        @warn "Conversion Failed: μ < tolerance."
        return fill(NaN, 6)
    end

    # Extract position and velocity vectors from Cartesian state
    r̄ = cart[1:3]
    v̄ = cart[4:6]
    twopi = 2π

    # Compute position and velocity magnitudes and test for singularity
    r = norm(r̄)
    v = norm(v̄)
    if r < tol || v < tol
        @warn "Conversion failed: Orbit is singular due to degenerate position or velocity vector."
        return fill(NaN, 6)
    end

    # Compute angular momentum and energy
    energy = v^2 / 2 - μ / r
    h̄ = cross(r̄, v̄)
    h = norm(h̄)
    if h < tol
        @warn "Conversion Failed: Orbit is singular due to degenerate angular momentum."
        return fill(NaN, 6)
    end

    # Compute inclination
    i = acos(clamp(h̄[3] / h, -1.0, 1.0))

    # Compute node vector
    z̄ = [0.0, 0.0, 1.0]  # Unit vector along Z-axis (inertial reference)
    n̄ = cross(z̄, h̄)
    n = norm(n̄)

    # Compute eccentricity vector and magnitude
    rdotv = dot(r̄, v̄)
    ē = ((v^2 - μ / r) * r̄ - rdotv * v̄) / μ
    e = norm(ē)

    # Semi-major axis (a) and semi-latus rectum (p)
    if abs(1 - e) > tol
        a = -μ / (2 * energy)  # For elliptical and hyperbolic orbits
        p = a * (1 - e^2)
    else
        a = Inf  # Parabolic case: semi-major axis undefined
        p = h^2 / μ
    end

    # Check for special cases and compute Ω, ω, and ν accordingly
    if i <= tol && e > tol
        # Elliptical equatorial orbit, set Ω = 0.0, measure ω from x-axis
        Ω = 0.0
        ω = acos(clamp(ē[1] / e, -1.0, 1.0))
        ω = ē[2] < 0 ? twopi - ω : ω
        ν = acos(clamp(dot(ē, r̄) / (e * r), -1.0, 1.0))
        ν = rdotv < 0 ? twopi - ν : ν
    elseif i > tol && e <= tol
        # Circular inclined orbit, set ω = 0.0, measure ν from ascending node
        Ω = acos(clamp(n̄[1] / n, -1.0, 1.0))
        Ω = n̄[2] < 0 ? twopi - Ω : Ω
        ω = 0.0
        ν = acos(clamp(dot(n̄ / n, r̄ / r), -1.0, 1.0))
        ν = r̄[3] < 0 ? twopi - ν : ν
    elseif i <= tol && e <= tol
        # Circular equatorial orbit, Ω = ω = 0.0
        Ω = 0.0
        ω = 0.0
        ν = acos(clamp(r̄[1] / r, -1.0, 1.0))
        ν = r̄[2] < 0 ? twopi - ν : ν
    else
        # General orbit (non-circular, inclined)
        Ω = acos(clamp(n̄[1] / n, -1.0, 1.0))
        Ω = n̄[2] < 0 ? twopi - Ω : Ω
        ω = acos(clamp(dot(n̄ / n, ē / e), -1.0, 1.0))
        ω = ē[3] < 0 ? twopi - ω : ω
        ν = acos(clamp(dot(ē / e, r̄ / r), -1.0, 1.0))
        ν = rdotv < 0 ? twopi - ν : ν
    end

    return [a, e, i, Ω, ω, ν]
end