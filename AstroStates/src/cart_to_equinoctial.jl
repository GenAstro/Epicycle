# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using LinearAlgebra

# Equinoctial elements are singular at i = π, where p and q grow as 1/(1 + cos i). Near it the
# conversion does not fail, it loses precision: about a kilometre of round-trip error at
# i = π - 1e-6. Orbits with 1 + cos i below this are refused rather than returned degraded; the
# margin it leaves is i > π - 1.4e-4 rad (0.008°), where the retrograde MEE set applies.
const _EQUINOCTIAL_RETROGRADE_TOL = 1e-8

"""
    cart_to_equinoctial(cart::AbstractVector{<:Real}, μ::Real; tol::Real = 1e-12)

Convert a Cartesian state vector to equinoctial orbital elements.

# Arguments
- `cart`: Cartesian state `[x, y, z, vx, vy, vz]`
- `μ::Real`: gravitational parameter [length³/time²]
- `tol::Real`: tolerance for singularity detection (default = 1e-12)

# Returns
Equinoctial state `[a, h, k, p, q, λ]`:
- `a` : semi-major axis [length], > 0
- `h` : e⋅sin(ω + Ω), the eccentricity vector along ĝ
- `k` : e⋅cos(ω + Ω), the eccentricity vector along f̂
- `p` : tan(i/2)⋅sin(Ω)
- `q` : tan(i/2)⋅cos(Ω)
- `λ` : mean longitude Ω + ω + M [rad], in [0, 2π)

# Notes
- Equinoctial elements here are defined for elliptic orbits only. A parabolic or hyperbolic
  orbit, a singular state, or an orbit within 0.008° of i = π (where the elements lose
  precision) logs a warning and returns `NaN`s.
- All angles are in radians. Units consistent with `μ`.
- Note that in most cases in states, h is the magnitude of angular momentum.  But, not
  for equinoctial elements.

# Examples
```julia
cart = [6778.0, 0.0, 0.0, 0.0, 7.66, 0.0]
equinoctial = cart_to_equinoctial(cart, 398600.4418)
```
"""
function cart_to_equinoctial(cart::AbstractVector{<:Real}, μ::Real; tol::Real = 1e-12)
    if length(cart) != 6
        error("Input must be a 6-element Cartesian state vector [x, y, z, vx, vy, vz] .")
    end
    T = float(promote_type(eltype(cart), typeof(μ)))

    r̄ = SVector{3,T}(cart[1], cart[2], cart[3])
    v̄ = SVector{3,T}(cart[4], cart[5], cart[6])

    r = norm(r̄)
    v = norm(v̄)

    if r < tol
        @warn "Conversion failed: Position magnitude r = $r less than tol."
        return fill(T(NaN), 6)
    end
    if μ < tol
        @warn "Conversion failed: Gravitational parameter μ = $μ less than tol."
        return fill(T(NaN), 6)
    end

    # Angular momentum vector and check for radial/degenerate orbit
    ang_mom_vec = cross(r̄, v̄)
    h_mag = norm(ang_mom_vec)
    if h_mag < tol
        @warn "Conversion failed: Angular momentum near zero (radial or degenerate orbit)."
        return fill(T(NaN), 6)
    end

    # Eccentricity vector
    ē = ((v^2 - μ / r) * r̄ - dot(r̄, v̄) * v̄) / μ

    # Specific energy and semi-major axis
    ξ = v^2 / 2 - μ / r
    a = -μ / (2 * ξ)
    if a < tol
        @warn "Conversion failed: Orbit is parabolic or hyperbolic (a = $a)."
        return fill(T(NaN), 6)
    end

    # Angular momentum unit vector; the elements are singular at i = π
    ĥ = ang_mom_vec / h_mag
    denom = 1 + ĥ[3]
    if denom < _EQUINOCTIAL_RETROGRADE_TOL
        @warn "Conversion failed: Equinoctial elements not defined for i ≈ π."
        return fill(T(NaN), 6)
    end

    # Equinoctial reference vectors f and g
    f̂ = SVector{3,T}(1 - ĥ[1]^2 / denom, -ĥ[1] * ĥ[2] / denom, -ĥ[1])
    ĝ = cross(ĥ, f̂)

    # Project eccentricity vector onto equinoctial frame
    h = dot(ē, ĝ)
    k = dot(ē, f̂)
    p = ĥ[1] / denom
    q = -ĥ[2] / denom

    # Eccentric longitude, then the mean longitude
    X1 = dot(r̄, f̂)
    Y1 = dot(r̄, ĝ)
    sqrt1 = sqrt(1 - h^2 - k^2)
    β = 1 / (1 + sqrt1)

    cosF = k + ((1 - k^2 * β) * X1 - h * k * β * Y1) / (a * sqrt1)
    sinF = h + ((1 - h^2 * β) * Y1 - h * k * β * X1) / (a * sqrt1)
    F = atan(sinF, cosF)

    λ = mod(F + h * cosF - k * sinF, 2 * T(π))

    return T[a, h, k, p, q, λ]
end
