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
    modkep_to_kep(modkep::Vector{<:Real}; tol::Float64 = 1e-12)

Convert a Modified Keplerian state to a classical Keplerian state.

# Arguments
- `modkep::Vector{<:Real}`: Modified Keplerian state vector `[rₚ, rₐ, i, Ω, ω, ν]`
    - `rₚ` : radius of periapsis
    - `rₐ` : radius of apoapsis 
    - `i`  : inclination 
    - `Ω`  : right ascension of ascending node 
    - `ω`  : argument of periapsis 
    - `ν`  : true anomaly

- `tol::Float64`: tolerance for singularity and consistency checks (default = 1e-12)

# Returns
A classical Keplerian state `[a, e, i, Ω, ω, ν]` or `fill(NaN, 6)` if invalid.

# Notes
- Returns `NaN` vector if input describes a parabolic, undefined, or inconsistent orbit.
- Assumes input angles are in radians and distances use consistent units.
"""
function modkep_to_kep(modkep::Vector{<:Real}; tol::Float64 = 1e-12)
    if length(modkep) != 6
        error("Input vector must contain six elements: [rₚ, rₐ, i, Ω, ω, ν]")
    end

    rₚ, rₐ, i, Ω, ω, ν = modkep

    # Check for singular orbits
    if abs(rₐ) < tol
        @warn "Conversion failed: Radius of apoapsis must not be zero."
        return fill(NaN, 6)
    end
    if rₐ < rₚ && rₐ > 0.0
        @warn "Conversion failed: Inconsistent Modified Keplerian state. If rₐ < rₚ, then rₐ must be negative."
        return fill(NaN, 6)
    end
    if rₚ <= tol || abs(rₐ) <= tol
        @warn "Conversion failed: Singular conic section."
        return fill(NaN, 6)
    end

    # Compute eccentricity and semi-major axis
    rp_by_ra = rₚ / rₐ
    e = (1 - rp_by_ra) / (1 + rp_by_ra)
    a = rₚ / (1 - e)

    return [a, e, i, Ω, ω, ν]
end
