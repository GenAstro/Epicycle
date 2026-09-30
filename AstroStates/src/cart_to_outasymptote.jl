# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using LinearAlgebra

"""
    cart_to_outasymptote(cart::AbstractVector{<:Real}, μ::Real; tol::Real = 1e-12) -> Vector

Convert a Cartesian state to outgoing asymptote parameters.

# Arguments
- `cart`: Cartesian state `[x, y, z, vx, vy, vz]`
- `μ::Real`: gravitational parameter
- `tol::Real`: tolerance for detecting singularities

# Returns
A 6-element vector: `[rₚ, C₃, λₐ, δₐ, θᵦ, ν]`
- rₚ : radius of periapsis
- C₃ : characteristic energy
- λₐ : right ascension of outgoing asymptote, in [0, 2π)
- δₐ : declination of outgoing asymptote, in [-π/2, π/2]
- θᵦ : B-plane angle, in [0, 2π)
- ν : true anomaly, in [0, 2π)

For an elliptic orbit (C₃ < 0) there is no asymptote, and the apoapsis direction is used in its
place. A circular or parabolic orbit, or an asymptote along the z-axis, logs a warning and
returns `NaN`s.

# Notes
- Angles are in radians.
- Dimensional quantities are consistent with units of μ.

# Examples
```julia
cart = [10000.0, 0.0, 0.0, 0.0, 12.0, 0.0]  # Hyperbolic trajectory
outasym = cart_to_outasymptote(cart, 398600.4418)
```
"""
cart_to_outasymptote(cart::AbstractVector{<:Real}, μ::Real; tol::Real = 1e-12) =
    _cart_to_asymptote(cart, μ, 1, tol)

# The outgoing (dir = 1) and incoming (dir = -1) asymptotes differ only in the sign of the
# in-plane component of the asymptote direction.
function _cart_to_asymptote(cart::AbstractVector{<:Real}, μ::Real, dir::Int, tol::Real)
    if length(cart) != 6
        error("Input must be a 6-element Cartesian vector: [x, y, z, vx, vy, vz]")
    end
    T = float(promote_type(eltype(cart), typeof(μ)))

    r̄ = SVector{3,T}(cart[1], cart[2], cart[3])
    v̄ = SVector{3,T}(cart[4], cart[5], cart[6])

    # Position, velocity and angular momentum magnitudes
    r = norm(r̄)
    v = norm(v̄)
    h̄ = cross(r̄, v̄)
    h = norm(h̄)

    # Degenerate cases
    if r < tol || v < tol
        @warn "Conversion failed: Orbit is singular due to degenerate position or velocity vector."
        return fill(T(NaN), 6)
    end
    if h < tol
        @warn "Conversion failed: Orbit is singular with zero angular momentum."
        return fill(T(NaN), 6)
    end

    # Eccentricity vector and magnitude
    ē = cross(v̄, h̄) / μ - r̄ / r
    e = norm(ē)
    if e < tol
        @warn "Conversion failed: Orbit is circular."
        return fill(T(NaN), 6)
    end

    # Characteristic energy; a parabola has no asymptote parameters
    C₃ = v^2 - 2μ / r
    if isapprox(C₃, 0; atol=tol)
        @warn "Conversion failed: Orbit is parabolic."
        return fill(T(NaN), 6)
    end

    # Radius of periapsis. rₚ = h²/(μ(1+e)) ≥ tol²/μ given the h check above, so it needs no test.
    a = -μ / C₃
    rₚ = a * (1 - e)

    # Asymptote unit vector ŝ; for an ellipse, the apoapsis direction
    if C₃ > tol
        fac = 1 / (1 + (C₃ * h^2) / μ^2)
        ŝ = fac * (dir * sqrt(C₃) / μ * cross(h̄, ē) - ē)
    else
        ŝ = -ē / e
    end

    # The B-plane axes need an asymptote off the z-axis
    sxy = hypot(ŝ[1], ŝ[2])
    if sxy < tol
        @warn "Conversion failed: Asymptote vector is aligned with the z-axis."
        return fill(T(NaN), 6)
    end

    # B-plane axes and the B-plane angle
    Ê = SVector{3,T}(-ŝ[2], ŝ[1], 0) / sxy          # ẑ × ŝ, normalised
    N̂ = cross(ŝ, Ê)
    b̄ = cross(h̄, ŝ)
    θᵦ = _wrap_2pi(atan(dot(b̄, Ê), -dot(b̄, N̂)))

    # Asymptote orientation angles
    δₐ = atan(ŝ[3], sxy)
    λₐ = _wrap_2pi(atan(ŝ[2], ŝ[1]))

    # True anomaly, from periapsis in the direction of motion
    ê = ē / e
    ν = _plane_angle(r̄, ê, cross(h̄ / h, ê))

    return T[rₚ, C₃, λₐ, δₐ, θᵦ, ν]
end
