# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

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

# Examples
```julia
modkep = [6778.0, 42164.0, π/6, 0.0, 0.0, 0.0]
kep = modkep_to_kep(modkep)
```
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
