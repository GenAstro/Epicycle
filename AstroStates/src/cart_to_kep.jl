# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    cart_to_kep(cart::AbstractVector{<:Real}, μ::Real; tol::Real=1e-12)

Convert a Cartesian state vector to Keplerian orbital elements.

# Arguments
- `cart`: Cartesian state `[x, y, z, vx, vy, vz]`
- `μ::Real`: Gravitational parameter
- `tol::Real`: Tolerance for detecting circular (`e ≤ tol`) and equatorial (`sin i ≤ tol`) orbits
  (default: `1e-12`)

# Returns
A vector `[a, e, i, Ω, ω, ν]` where:
- `a`: semi-major axis; `Inf` for a parabolic orbit
- `e`: eccentricity
- `i`: inclination, in [0, π]
- `Ω`: right ascension of ascending node (RAAN), in [0, 2π)
- `ω`: argument of periapsis, in [0, 2π)
- `ν`: true anomaly, in [0, 2π)

The element type follows the inputs, so `Float32` and `ForwardDiff.Dual` states stay so. A
singular state (zero position, velocity or angular momentum, or μ below `tol`) logs a warning and
returns `NaN`s.

# Notes
- Angles are in radians; dimensional quantities must use units consistent with μ.
- Where an element is undefined it is pinned, as GMAT does:
  * equatorial orbits (i = 0 or π): Ω = 0, and ω is the longitude of periapsis, measured from +x
    in the direction of motion;
  * circular orbits: ω = 0, and ν is the argument of latitude, measured from the ascending node;
  * circular equatorial orbits: Ω = ω = 0, and ν is the true longitude, measured from +x in the
    direction of motion.
- Angles are recovered with `atan(y, x)`, so they keep full precision near 0 and π and have
  derivatives there.

# Examples
```julia
cart = [6778.0, 0.0, 0.0, 0.0, 7.66, 0.0]
kep = cart_to_kep(cart, 398600.4418)
```
"""
function cart_to_kep(cart::AbstractVector{<:Real}, μ::Real; tol::Real=1e-12)
    if length(cart) != 6
        error("Input vector must have exactly six elements: [x, y, z, vx, vy, vz].")
    end
    T = float(promote_type(eltype(cart), typeof(μ)))

    if μ < tol
        @warn "Conversion Failed: μ < tolerance."
        return fill(T(NaN), 6)
    end

    r̄ = SVector{3,T}(cart[1], cart[2], cart[3])
    v̄ = SVector{3,T}(cart[4], cart[5], cart[6])

    # Position and velocity magnitudes, and the singular states
    r = norm(r̄)
    v = norm(v̄)
    if r < tol || v < tol
        @warn "Conversion failed: Orbit is singular due to degenerate position or velocity vector."
        return fill(T(NaN), 6)
    end

    # Angular momentum and energy
    energy = v^2 / 2 - μ / r
    h̄ = cross(r̄, v̄)
    h = norm(h̄)
    if h < tol
        @warn "Conversion Failed: Orbit is singular due to degenerate angular momentum."
        return fill(T(NaN), 6)
    end

    # Eccentricity vector and magnitude
    ē = ((v^2 - μ / r) * r̄ - dot(r̄, v̄) * v̄) / μ
    e = norm(ē)

    # Semi-major axis; undefined (infinite) for a parabola
    a = abs(1 - e) > tol ? -μ / (2 * energy) : T(Inf)

    i, Ω, ω, p̂, q̂ = _orbit_orientation(h̄, ē, e, tol)
    ν = _plane_angle(r̄, p̂, q̂)

    return T[a, e, i, Ω, ω, ν]
end
