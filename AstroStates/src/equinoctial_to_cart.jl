# Copyright (C) 2025 Gen Astro LLC
SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

using LinearAlgebra

"""
    equinoctial_to_cart(eq::Vector{<:Real}, μ::Real; tol::Float64 = 1e-12)

Convert equinoctial elements to Cartesian state.

# Arguments
- `eq::Vector{<:Real}`: Equinoctial state `[a, h, k, p, q, λ]`
    - `a` : semi-major axis [length]
    - `h` : e⋅g component of eccentricity vector
    - `k` : e⋅f component of eccentricity vector
    - `p` : tan(i/2)⋅cos(Ω)
    - `q` : tan(i/2)⋅sin(Ω)
    - `λ` : mean longitude [rad]
- `μ::Real`: gravitational parameter [length³/time²]
- `tol::Float64`: numerical tolerance (default: 1e-12)

# Returns
Cartesian state `[x, y, z, vx, vy, vz]`

# Notes
- Returns `fill(NaN, 6)` if eccentricity is too high or radius becomes non-physical.
- Assumes all angles are in radians and other units are consistent with μ.

# Examples
```julia
eq = [7000.0, 0.01, 0.0, 0.1, 0.0, π/4]
cart = equinoctial_to_cart(eq, 398600.4418)
```
"""
function equinoctial_to_cart(eq::Vector{<:Real}, μ::Real; tol::Float64 = 1e-12)
    if length(eq) != 6
        error("Input vector must contain six equinoctial elements: [a, h, k, p, q, λ]")
    end

    a, h, k, p, q, λ = eq
    e = sqrt(h^2 + k^2)

    if e > 1 - tol
        @warn "Conversion failed: Eccentricity (e = $e) exceeds bound for equinoctial formulation."
        return fill(NaN, 6)
    end

    # Solve for true longitude F from mean longitude λ using Newton-Raphson
    F = λ
    converged = false
    for i in 1:100
        fF = F + h * cos(F) - k * sin(F) - λ
        f′F = 1 - h * sin(F) - k * cos(F)
        ΔF = -fF / f′F
        F += ΔF
        if abs(ΔF) < tol
            converged = true
            break
        end
    end

    if !converged
        @warn "Conversion failed: Reached max iterations in conversion from mean to eccentric longitude."
        return fill(NaN, 6)
    end

    F = F < 0 ? F + 2π : F

    # Compute r, beta
    sqrt1 = sqrt(1 - h^2 - k^2)
    β = 1 / (1 + sqrt1)
    n = sqrt(μ / a^3)
    cosF, sinF = cos(F), sin(F)
    r = a * (1 - k * cosF - h * sinF)

    if r <= 0
        @warn "Conversion failed: Radius (r = $r) is non-physical."
        return fill(NaN, 6)
    end

    # Position and velocity in orbital plane
    X₁ = a * ((1 - h^2 * β) * cosF + h * k * β * sinF - k)
    Y₁ = a * ((1 - k^2 * β) * sinF + h * k * β * cosF - h)

    Ẋ₁ = (n * a^2 / r) * (h * k * β * cosF - (1 - h^2 * β) * sinF)
    Ẏ₁ = (n * a^2 / r) * ((1 - k^2 * β) * cosF - h * k * β * sinF)

    # Construct Q matrix and extract f̂, ĝ directions
    Q = [
        1 - p^2 + q^2      2p*q               2p;
        2p*q               1 + p^2 - q^2     -2q;
       -2p                 2q                 1 - p^2 - q^2
    ]

    Q2 = Q / (1 + p^2 + q^2)

    f̂ = normalize(Q2[:, 1])
    ĝ = normalize(Q2[:, 2])

    # Final position and velocity vectors
    r̄ = X₁ * f̂ + Y₁ * ĝ
    v̄ = Ẋ₁ * f̂ + Ẏ₁ * ĝ

    return [r̄; v̄]
end
