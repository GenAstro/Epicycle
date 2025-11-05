# Copyright (C) 2025 Gen Astro LLC
SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    alt_equinoctial_to_equinoctial(alt::Vector{<:Real}; tol::Float64 = 1e-12)

Convert alternate equinoctial elements to standard equinoctial elements.

# Arguments
- `alt::Vector{<:Real}`: Alternate equinoctial state `[a, h, k, altp, altq, λ]`
    - `a` : semi-major axis [length]
    - `h` : e⋅g component of eccentricity vector
    - `k` : e⋅f component of eccentricity vector
    - `altp` : sin(i/2)⋅cos(Ω)
    - `altq` : sin(i/2)⋅sin(Ω)
    - `λ` : mean longitude [rad]
- `tol::Float64`: tolerance for inclination singularity (default = 1e-12)

# Returns
Standard equinoctial state `[a, h, k, p, q, λ]`
    - `a`  : semi-major axis [length]
    - `h`  : e⋅g component of eccentricity vector
    - `k`  : e⋅f component of eccentricity vector
    - `p`  : tan(i/2)⋅cos(Ω)
    - `q`  : tan(i/2)⋅sin(Ω)
    - `λ`  : mean longitude [rad]

# Notes
- Fails if inclination approaches 180°.
- All angles in radians. Units consistent with `μ`.

# Examples
```julia
alt_eq = [7000.0, 0.01, 0.0, 0.05, 0.0, π/4]
std_eq = alt_equinoctial_to_equinoctial(alt_eq)
```
"""
function alt_equinoctial_to_equinoctial(alt::Vector{<:Real}; tol::Float64 = 1e-12)
    if length(alt) != 6
        error("Input must be a 6-element alternate equinoctial vector: [a, h, k, altp, altq, λ]")
    end

    a, h, k, altp, altq, λ = alt

    # Compute inclination from alternate p/q
    i = 2 * asin(clamp(sqrt(altp^2 + altq^2),-1.0,1.0))

    if abs(i - π) < tol
        @warn "Conversion failed: Inclination i ≈ π (180°) causes singularity in equinoctial form."
        return fill(NaN, 6)
    end

    p = altp / cos(i / 2)
    q = altq / cos(i / 2)

    return [a, h, k, p, q, λ]
end
