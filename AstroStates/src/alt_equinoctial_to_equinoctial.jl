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
