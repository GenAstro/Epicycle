# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    alt_equinoctial_to_equinoctial(alt::AbstractVector{<:Real};
        tol::Real = _ALT_EQUINOCTIAL_RETROGRADE_TOL)

Convert alternate equinoctial elements to standard equinoctial elements.

# Arguments
- `alt::AbstractVector{<:Real}`: Alternate equinoctial state `[a, h, k, altp, altq, λ]`
    - `a` : semi-major axis [length]
    - `h` : e⋅g component of eccentricity vector
    - `k` : e⋅f component of eccentricity vector
    - `altp` : sin(i/2)⋅sin(Ω)
    - `altq` : sin(i/2)⋅cos(Ω)
    - `λ` : mean longitude Ω + ω + M [rad]
- `tol::Real`: how close to i = π, in radians, the elements are refused (default 1.4e-4 rad,
  about 0.008°, where the equinoctial set refuses too; nearer, the conversion loses precision)

# Returns
Standard equinoctial state `[a, h, k, p, q, λ]`
    - `a`  : semi-major axis [length]
    - `h`  : e⋅g component of eccentricity vector
    - `k`  : e⋅f component of eccentricity vector
    - `p`  : tan(i/2)⋅sin(Ω)
    - `q`  : tan(i/2)⋅cos(Ω)
    - `λ`  : mean longitude Ω + ω + M [rad]

# Notes
- Fails if inclination approaches 180°.
- All angles in radians. Units consistent with `μ`.

# Examples
```julia
alt_eq = [7000.0, 0.01, 0.0, 0.05, 0.0, π/4]
std_eq = alt_equinoctial_to_equinoctial(alt_eq)
```
"""
function alt_equinoctial_to_equinoctial(alt::AbstractVector{<:Real};
        tol::Real = _ALT_EQUINOCTIAL_RETROGRADE_TOL)
    if length(alt) != 6
        error("Input must be a 6-element alternate equinoctial vector: [a, h, k, altp, altq, λ]; got $(length(alt)) elements.")
    end
    _has_nan(alt) && return fill(float(eltype(alt))(NaN), 6)

    a, h, k, altp, altq, λ = alt

    # altp² + altq² = sin²(i/2), so cos²(i/2) = 1 - altp² - altq². Written that way there is no
    # sqrt(altp² + altq²), which has no derivative at an equatorial orbit, where the elements do.
    # i within tol of π is cos²(i/2) below (tol / 2)².
    cos2_half_i = 1 - altp^2 - altq^2
    # A value a rounding error below zero is an orbit at i = π, refused just below; only one
    # clearly below zero has no inclination
    if cos2_half_i < -4 * eps(float(eltype(alt)))
        @warn "Conversion failed: altp² + altq² = $(altp^2 + altq^2) exceeds 1; it is sin²(i/2), so no inclination has it."
        return fill(float(eltype(alt))(NaN), 6)
    end
    if cos2_half_i < (tol / 2)^2
        @warn "Conversion failed: Inclination i ≈ π (180°) causes singularity in equinoctial form."
        return fill(float(eltype(alt))(NaN), 6)
    end

    cos_half_i = sqrt(cos2_half_i)
    p = altp / cos_half_i
    q = altq / cos_half_i

    return [a, h, k, p, q, λ]
end
