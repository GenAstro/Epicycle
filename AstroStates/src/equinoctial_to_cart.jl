# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using LinearAlgebra

"""
    equinoctial_to_cart(eq::AbstractVector{<:Real}, μ::Real; tol::Real = 1e-12)

Convert equinoctial elements to Cartesian state.

# Arguments
- `eq`: Equinoctial state `[a, h, k, p, q, λ]`
    - `a` : semi-major axis [length], > 0
    - `h` : e⋅sin(ω + Ω)
    - `k` : e⋅cos(ω + Ω)
    - `p` : tan(i/2)⋅sin(Ω)
    - `q` : tan(i/2)⋅cos(Ω)
    - `λ` : mean longitude Ω + ω + M [rad]
- `μ::Real`: gravitational parameter [length³/time²]
- `tol::Real`: numerical tolerance (default: 1e-12)

# Returns
Cartesian state `[x, y, z, vx, vy, vz]`

# Notes
- Equinoctial elements here are defined for elliptic orbits only. `a ≤ 0`, `e ≥ 1`, or a
  non-positive μ logs a warning and returns `fill(NaN, 6)`.
- Kepler's equation is solved by a safeguarded Newton iteration, which converges for any
  eccentricity below 1.
- Assumes all angles are in radians and other units are consistent with μ.

# Examples
```julia
eq = [7000.0, 0.01, 0.0, 0.1, 0.0, π/4]
cart = equinoctial_to_cart(eq, 398600.4418)
```
"""
function equinoctial_to_cart(eq::AbstractVector{<:Real}, μ::Real; tol::Real = 1e-12)
    if length(eq) != 6
        error("Input vector must contain six equinoctial elements: [a, h, k, p, q, λ]")
    end
    T = float(promote_type(eltype(eq), typeof(μ)))

    a, h, k, p, q, λ = eq
    e = sqrt(h^2 + k^2)

    if μ < tol
        @warn "Conversion failed: Gravitational parameter μ = $μ less than tol."
        return fill(T(NaN), 6)
    end
    if a <= tol
        @warn "Conversion failed: Semi-major axis a = $a; equinoctial elements are defined for " *
              "elliptic orbits only."
        return fill(T(NaN), 6)
    end
    if e > 1 - tol
        @warn "Conversion failed: Eccentricity (e = $e) exceeds bound for equinoctial formulation."
        return fill(T(NaN), 6)
    end

    # Eccentric longitude F from the mean longitude. With F = φ + E, φ = atan(h, k) the longitude
    # of periapsis, λ = F + h cos F - k sin F is Kepler's equation λ - φ = E - e sin E.
    φ = atan(h, k)
    F = φ + _eccentric_anomaly(λ - φ, e)

    # Radius
    sqrt1 = sqrt(1 - h^2 - k^2)
    β = 1 / (1 + sqrt1)
    n = sqrt(μ / a^3)
    sinF, cosF = sincos(F)
    r = a * (1 - k * cosF - h * sinF)

    # Position and velocity in orbital plane
    X₁ = a * ((1 - h^2 * β) * cosF + h * k * β * sinF - k)
    Y₁ = a * ((1 - k^2 * β) * sinF + h * k * β * cosF - h)

    Ẋ₁ = (n * a^2 / r) * (h * k * β * cosF - (1 - h^2 * β) * sinF)
    Ẏ₁ = (n * a^2 / r) * ((1 - k^2 * β) * cosF - h * k * β * sinF)

    # Equinoctial frame: the first two columns of Q / (1 + p² + q²), which are unit vectors
    s = 1 + p^2 + q^2
    f̂ = SVector{3,T}(1 - p^2 + q^2, 2p * q, -2p) / s
    ĝ = SVector{3,T}(2p * q, 1 + p^2 - q^2, 2q) / s

    # Final position and velocity vectors
    r̄ = X₁ * f̂ + Y₁ * ĝ
    v̄ = Ẋ₁ * f̂ + Ẏ₁ * ĝ

    return T[r̄[1], r̄[2], r̄[3], v̄[1], v̄[2], v̄[3]]
end
