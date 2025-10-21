# Copyright 2025 Gen Astro LLC. All Rights Reserved.
#
# This software is licensed under the GNU AGPL v3.0,
# WITHOUT ANY WARRANTY, including implied warranties of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# This file may also be used under a commercial license,
# if one has been purchased from Gen Astro LLC.
#
# By modifying this software, you agree to the terms of the
# Gen Astro LLC Contributor License Agreement.

"""
    equinoctial_to_alt_equinoctial(eq::Vector{<:Real}; tol::Float64 = 1e-12)

Convert equinoctial elements to alternate equinoctial elements.

# Arguments
- `eq::Vector{<:Real}`: Equinoctial state `[a, h, k, p, q, λ]`
    - `a`  : semi-major axis [length]
    - `h`  : e⋅g component of eccentricity vector
    - `k`  : e⋅f component of eccentricity vector
    - `p`  : tan(i/2)⋅cos(Ω)
    - `q`  : tan(i/2)⋅sin(Ω)
    - `λ`  : mean longitude [rad]
- `tol::Float64`: tolerance for inclination singularity detection (default: 1e-12)

# Returns
Alternate equinoctial state `[a, h, k, altp, altq, λ]`
    - `a` : semi-major axis [length]
    - `h` : e⋅g component of eccentricity vector
    - `k` : e⋅f component of eccentricity vector
    - `altp` : sin(i/2)⋅cos(Ω)
    - `altq` : sin(i/2)⋅sin(Ω)
    - `λ` : mean longitude [rad]

# Notes
- Singular when inclination is near 180°, returns `fill(NaN, 6)` in that case.
- Units and angles must be consistent; λ in radians.
"""
function equinoctial_to_alt_equinoctial(eq::Vector{<:Real}; tol::Float64 = 1e-12)
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

    altp = p * cos(i / 2)
    altq = q * cos(i / 2)

    return [a, h, k, altp, altq, λ]
end
