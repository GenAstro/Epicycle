# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    equinoctial_to_alt_equinoctial(eq::AbstractVector{<:Real}; tol::Real = 1e-12)

Convert equinoctial elements to alternate equinoctial elements.

# Arguments
- `eq::AbstractVector{<:Real}`: Equinoctial state `[a, h, k, p, q, λ]`
    - `a`  : semi-major axis [length]
    - `h`  : e⋅g component of eccentricity vector
    - `k`  : e⋅f component of eccentricity vector
    - `p`  : tan(i/2)⋅sin(Ω)
    - `q`  : tan(i/2)⋅cos(Ω)
    - `λ`  : mean longitude Ω + ω + M [rad]
- `tol::Real`: tolerance for inclination singularity detection (default: 1e-12)

# Returns
Alternate equinoctial state `[a, h, k, altp, altq, λ]`
    - `a` : semi-major axis [length]
    - `h` : e⋅g component of eccentricity vector
    - `k` : e⋅f component of eccentricity vector
    - `altp` : sin(i/2)⋅sin(Ω)
    - `altq` : sin(i/2)⋅cos(Ω)
    - `λ` : mean longitude Ω + ω + M [rad]

# Notes
- Singular when inclination is near 180°, returns `fill(NaN, 6)` in that case.
- Units and angles must be consistent; λ in radians.

# Examples
```julia
eq = [7000.0, 0.01, 0.0, 0.1, 0.0, π/4]
alt_eq = equinoctial_to_alt_equinoctial(eq)
```
"""
function equinoctial_to_alt_equinoctial(eq::AbstractVector{<:Real}; tol::Real = 1e-12)
    if length(eq) != 6
        error("Input must be a 6-element equinoctial vector: [a, h, k, p, q, λ]")
    end

    a, h, k, p, q, λ = eq

    # Compute inclination from equinoctial p, q
    i = 2 * atan(sqrt(p^2 + q^2))

    if abs(i - π) < tol
        @warn "Conversion failed: Inclination i ≈ π (180°) causes singularity in alternate equinoctial form."
        return fill(NaN, 6)
    end

    # cos(i/2) from p and q directly: through i it has no derivative at i = 0
    cos_half_i = 1 / sqrt(1 + p^2 + q^2)
    altp = p * cos_half_i
    altq = q * cos_half_i

    return [a, h, k, altp, altq, λ]
end
