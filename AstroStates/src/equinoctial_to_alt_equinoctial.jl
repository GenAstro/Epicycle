# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    equinoctial_to_alt_equinoctial(eq::AbstractVector{<:Real};
        tol::Real = _ALT_EQUINOCTIAL_RETROGRADE_TOL)

Convert equinoctial elements to alternate equinoctial elements.

# Arguments
- `eq::AbstractVector{<:Real}`: Equinoctial state `[a, h, k, p, q, λ]`
    - `a`  : semi-major axis [length]
    - `h`  : e⋅g component of eccentricity vector
    - `k`  : e⋅f component of eccentricity vector
    - `p`  : tan(i/2)⋅sin(Ω)
    - `q`  : tan(i/2)⋅cos(Ω)
    - `λ`  : mean longitude Ω + ω + M [rad]
- `tol::Real`: how close to i = π, in radians, the elements are refused (default 1.4e-4 rad,
  about 0.008°, where the equinoctial set refuses too; nearer, the conversion loses precision)

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
function equinoctial_to_alt_equinoctial(eq::AbstractVector{<:Real};
        tol::Real = _ALT_EQUINOCTIAL_RETROGRADE_TOL)
    if length(eq) != 6
        error("Input must be a 6-element equinoctial vector: [a, h, k, p, q, λ]; got $(length(eq)) elements.")
    end
    _has_nan(eq) && return fill(float(eltype(eq))(NaN), 6)

    a, h, k, p, q, λ = eq

    # p² + q² = tan²(i/2), so cos(i/2) = 1 / sqrt(1 + p² + q²). Written that way there is no
    # sqrt(p² + q²), which has no derivative at an equatorial orbit, where the elements do.
    # i within tol of π is p² + q² above (2 / tol)².
    t2 = p^2 + q^2
    if t2 > (2 / tol)^2
        @warn "Conversion failed: Inclination i ≈ π (180°) causes singularity in alternate equinoctial form."
        return fill(float(eltype(eq))(NaN), 6)
    end

    cos_half_i = 1 / sqrt(1 + t2)
    altp = p * cos_half_i
    altq = q * cos_half_i

    return [a, h, k, altp, altq, λ]
end
